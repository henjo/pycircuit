"""`PAC`, the periodic small-signal and noise analyses over a `PSS` operating
point, and `SidebandResponse`.
"""
import numpy as np
import warnings
from pycircuit.circuit.analysis import Analysis
from pycircuit.circuit.analysis import Parameter
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.circuit import gnd
import pycircuit.circuit.analysis as analysis
from ._numerics import _arnoldi_gmres
from ._numerics import freq_analysis
from ._numerics import periodic_spline_weights
from .events import EventColumns


class SidebandResponse(object):
    """Every input band that lands on ONE output frequency.

    ⚠ THE POINT OF THIS OBJECT IS THAT PAC'S ANSWER IS NOT ONE NUMBER.
    Kundert: *"for a single output frequency there may be many transfer
    functions from a single input"*.  Conversion gain is one entry; image
    rejection, LO feedthrough and supply rejection are others, and a
    caller who takes a single coefficient and calls it "the gain" has
    silently picked one of them.

    ⚠⚠ AND THE BANDS SIT AT DIFFERENT INPUT FREQUENCIES, WHICH IS THE
    THING THAT IS EASY TO GET WRONG.  `adjoint_sideband_row` is indexed by
    the INPUT frequency: a source at `f` reaches the output at
    `f + l f0` through sideband `l`.  Fixing the OUTPUT instead means each
    sideband is fed from its own input band, `f_in = f_out - l f0`.  So
    image rejection is NOT `H_l` against `H_-l` at one input frequency --
    it is two different input bands mapping onto one output, and computing
    it the first way gives a plausible number for a different quantity.

    Attributes: `f_out`, `f0`, `sidebands`, `inputs` (the `f_in` per
    sideband) and `rows` (`(len(sidebands), m)`, source-indexed).
    """

    __slots__ = ('f_out', 'f0', 'sidebands', 'inputs', 'rows')

    def __init__(self, f_out, f0, sidebands, inputs, rows):
        self.f_out = float(f_out)
        self.f0 = float(f0)
        self.sidebands = list(sidebands)
        self.inputs = list(inputs)
        self.rows = np.asarray(rows)

    def _index(self, l):
        try:
            return self.sidebands.index(int(l))
        except ValueError:
            raise KeyError(
                'sideband %r was not computed; this response carries %r'
                % (l, self.sidebands))

    def input_frequency(self, l):
        """The input band feeding sideband `l`: `f_out - l f0`."""
        return self.inputs[self._index(l)]

    def transfer(self, l, source):
        """The complex coefficient from `source` at `f_in(l)` to `f_out`."""
        return complex(self.rows[self._index(l)][int(source)])

    def gain_db(self, l, source):
        """`20 log10 |H_l|` -- one band's conversion gain, named as such."""
        mag = abs(self.transfer(l, source))
        if mag == 0.0:
            return -np.inf
        return 20.0 * np.log10(mag)

    def rejection_db(self, wanted, other, source):
        """How far `other` sits below `wanted`, in dB.

        ⚠ WHICH REJECTION THIS IS DEPENDS ENTIRELY ON WHICH TWO SIDEBANDS
        ARE NAMED, and the method refuses to guess.  Image rejection is
        the wanted band against the one mirrored about the LO; LO
        feedthrough is the wanted band against `l` such that `f_in = 0`.
        Naming them at the call site is the difference between a number
        and a labelled number.
        """
        w = abs(self.transfer(wanted, source))
        o = abs(self.transfer(other, source))
        if o == 0.0:
            return np.inf
        if w == 0.0:
            return -np.inf
        return 20.0 * np.log10(w / o)

    def __repr__(self):
        return ('SidebandResponse(f_out=%g, f0=%g, sidebands=%r)'
                % (self.f_out, self.f0, self.sidebands))


def _output_weights(output, width):
    """The complex output functional of an adjoint row, `width` wide: an
    index is a unit vector, an array its own entries, zero-padded."""
    d = np.zeros(width, dtype=complex)
    if np.isscalar(output):
        d[int(output)] = 1.0
    else:
        out = np.asarray(output, dtype=complex).ravel()
        d[:len(out)] = out
    return d


def _output_row(output, m):
    """The real output row of a spectrum: an index (a 0-d value) is a unit
    vector of width `m`, an array its first `m` entries."""
    d = np.asarray(output)
    if d.ndim == 0:
        row = np.zeros(m, dtype=float)
        row[int(d)] = 1.0
        return row
    return np.asarray(d, dtype=float).ravel()[:m]


class PAC(Analysis):
    """Small-signal analysis over a periodic operating point, matrix-free.

    The operator is the monodromy.  The periodic small-signal system is

        (L + alpha B) v = -u,   alpha = exp(-2j pi f T)

    with `L` the block lower bidiagonal discretisation over the period and
    `B` the periodic wrap.  Telichevesky, Kundert & White (DAC 1996) use
    `L^-1` as a preconditioner:

        (I + alpha L^-1 B) v = -L^-1 u

    Applying `L^-1` is forward substitution through the timesteps -- the
    recursion PSS already runs against stored factors -- and `B` is
    confined to the first `m` rows and last `m` columns, so `L^-1 B` acts
    only on the LAST block
    (`test_the_pac_operator_is_the_monodromy_and_L_is_never_formed`).
    What is left is `m x m`:

        (I - alpha M) y_0 = alpha w(f)

    with `M` the monodromy and `w` the forced response over one period from
    a zero initial state.  `y_0` is the small-signal state at `t = 0`; one
    more driven replay gives the rest of the period.

    Neither `L` nor `B` is formed (`(N m)^2` complex entries, ~420 GiB at
    `N = 137`, `m = 1000`): the stored per-step factors are the
    preconditioner, and the only dense object is `m x m`, formed only on
    request.

    ⚠ Do not rebuild `L` with two terms per row: a two-step method's
    variational system has three, and a backward-Euler-shaped `L` for
    `trap` or `gear` is the operator of a different recursion (spectral
    radius 0 against the analytic 0.8546,
    `test_the_pac_L_is_backward_euler_only`).  `M` taken from the
    traversal cannot make that mistake: every step carries its own
    `(alphas, b)`.

    History: `doc/shooting_history.md`, `PAC`.
    """

    parameters  = [Parameter(name='analysis', desc='Analysis name',
                             default='ac')]

    ## How hard GMRES is asked to solve the `m x m` system, relative to the
    ## PSS reltol that produced the operating point.  Looser than the
    ## operating point itself would be answering a question the trajectory
    ## cannot support; much tighter buys nothing, because the linearisation
    ## is only as good as the trajectory.
    KRYLOV_FACTOR = 1e-2

    def __init__(self, cir, toolkit=None, **kvargs):
        self.parameters = super(PAC, self).parameters + self.parameters
        super(PAC, self).__init__(cir, toolkit=toolkit, **kvargs)

    def solve(self, pss, freqs, refnode=gnd, recycle=True):
        """Sideband response at each frequency in `freqs`.

        `pss` must be a CONVERGED `PSS` -- the periodic operating point is
        what this linearises about, and there is no meaningful small-signal
        answer about a non-solution.  `PSS.factored_period()` enforces it.

        `recycle` shares one Krylov subspace across the sweep, which is
        where the sweep's cost goes; see `_solve_subspace`.

        History: `doc/shooting_history.md`, `PAC.solve`.
        """
        toolkit = self.toolkit
        freqs = np.atleast_1d(np.asarray(freqs, dtype=float))
        ## the map on the state (a GLM's own, `_GLMPeriod.state_map`)
        fp = pss._state_map()
        T = float(fp.T)
        m = self.cir.n - 1

        irefnode = self.cir.get_node_index(refnode)
        if irefnode != pss.irefnode:
            raise ValueError(
                'PAC: refnode (index %d) differs from the PSS the operating '
                'point came from (index %d). The monodromy eliminated one '
                'and this would report against the other.'
                % (irefnode, pss.irefnode))
        ## ⚠ `analysis=` BY KEYWORD.  `Circuit.u(t, epar, analysis, ...)`
        ## takes `epar` second: positionally, 'ac' becomes the element
        ## parameter set and the TRANSIENT source vector comes back -- zero
        ## at `t = 0` for every sinusoid, so the analysis silently returns
        ## zeros.
        (u_ac,) = remove_row_col((self.cir.u(0, analysis=self.par.analysis),),
                                 irefnode, toolkit)
        if not np.any(np.asarray(u_ac)):
            raise ValueError(
                'PAC: the %r source vector is identically zero, so there is '
                'nothing to analyse. Independent sources take their '
                'small-signal amplitude from `vac`/`iac`, not from `va` -- '
                'a source with va= set and vac=0 drives the operating point '
                'and not this.' % self.par.analysis)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()

        ## the forced response at each frequency -- one period replay each,
        ## and unavoidable: the source is what changes across the sweep
        ## ⚠ THE MANUFACTURING STEP IS NOT IN `steps`, AND IT COSTS AN
        ## ORDER.  On the plain path the factored walk takes one
        ## step outside the loop to manufacture a history and folds it into
        ## the `opening` triple as a flat-history assumption.  The source is
        ## never applied at that step, so the driven response is first order
        ## in h whatever the method's order (trap: 2.00x per doubling plain,
        ## 4.00x with x0_unknown=True; euler unchanged).  Gear-2 takes the
        ## solved-history path and has no manufacturing step.
        if (fp.is_plain and not fp.open_at_x0
                and pss.par.method != 'euler'):
            warnings.warn(
                'PAC: this operating point was solved on the PLAIN path '
                'with a manufacturing step (method=%r, x0_unknown=False). '
                'The manufacturing step carries no small-signal source, so '
                'the response is FIRST order in the timestep whatever the '
                "method's own order -- measured 2.00x per doubling against "
                '4.00x for the same run with x0_unknown=True. The answer is '
                'not wrong, it is one order less accurate than the '
                'trajectory it came from. Re-solve with x0_unknown=True, or '
                "with method='gear', to get the method's own order."
                % pss.par.method,
                RuntimeWarning, stacklevel=2)

        self._check_circuit(pss)
        for f in freqs:
            self._check_harmonic(pss, f, 'a sweep point')

        ys, dthetas = self._forced_responses(pss, fp, freqs, u_ac, recycle)
        self.event_shifts = list(dthetas)
        outfreq, outV = [], []
        ## the complex time-domain response per frequency, `(times, y)` with
        ## `y` the response to the source `u_ac e^{j w t}` at the grid's
        ## nodes (reduced state); the sideband coefficients below are its
        ## periodic envelope's Fourier coefficients by the period quadrature,
        ## which on a strongly non-uniform grid is not an interpolating basis
        ## -- a time-domain reading comes from here, not from summing them
        self.time_response = []
        for f, y in zip(freqs, ys):
            self.time_response.append((np.asarray(fp.times, dtype=float)[:len(y)], y.copy()))
            ## `v(t) = y(t) exp(-j w t)` is T-periodic; its DFT is the
            ## sideband set
            tms = np.asarray(fp.times, dtype=float)[:len(y)]
            v = y * np.exp(-2j * np.pi * f * tms)[:, None]
            ## ⚠ `fp.times` spans `[0, T]` INCLUSIVE: the repeated endpoint is
            ## dropped before the DFT (kept, it puts the sidebands at
            ## `f0 (N-1)/N` and costs an order).  Guarded on the window, since
            ## the plain path's `[:len(y)]` need not be inclusive.
            ## ⚠ A NEGATIVE sideband frequency is folded to positive and its
            ## coefficient CONJUGATED -- that is the physical response there.
            ## Both are invisible on a circuit whose `v(t)` is constant over
            ## the period.
            if len(tms) > 1 and np.isclose(tms[-1] - tms[0], T,
                                           rtol=1e-9, atol=0.0):
                v, tms = v[:-1], tms[:-1]
            _wq = pss._period_quadrature(fp)
            if _wq is not None and len(_wq) == len(tms):
                ## non-uniform grid: the SAME trapezoid weights the adjoint
                ## inject uses, at the true times -- see `_period_quadrature`
                _ks = np.fft.fftshift(np.fft.fftfreq(len(tms), d=1.0 / len(tms)))
                sb = _ks / T
                V = np.tensordot(np.exp(-2j * np.pi * np.outer(
                    _ks, (tms - tms[0]) / T)) * _wq[None, :], v, axes=(1, 0))
            else:
                sb, V = freq_analysis(v, tms, axis=0)
            fs = np.asarray(sb, dtype=float) + f
            V = np.asarray(V)
            neg = fs < 0.0
            if np.any(neg):
                V = V.copy()
                V[neg] = np.conj(V[neg])
            outfreq.extend(np.abs(fs).tolist())
            outV.extend(V.tolist())

        order = np.argsort(np.asarray(outfreq))
        fout = np.asarray(outfreq)[order]
        X = np.asarray(outV)[order]
        X = np.concatenate((X[:, :irefnode],
                            np.zeros((len(fout), 1)),
                            X[:, irefnode:]), axis=1)
        self.result = analysis.CircuitResult(
            self.cir, x=X.T, xdot=None, sweep_values=fout,
            sweep_label='freq', sweep_unit='Hz')
        return self.result

    def _forced_responses(self, pss, fp, freqs, u_ac, recycle=True,
                          u_points=None):
        """The steady forced response to the source ``u_ac e^{j w t}`` (or a
        MODULATED one, `u_points`: `PSS._forced_replay`) at each frequency:
        per frequency the complex response at every node of the grid, `(N +
        1, m)`, and the crossings' shifts (None without state events).
        The operator solves (deflated on an oscillator, recycled across the
        sweep otherwise), the bordered event rows on a staged solve and the
        fixed-time correction of the node responses -- `PAC.solve`'s core
        (factored out 2026-09-25 so `_coloured_covariance` reads the same
        responses).  Sets `deflated` and `matvecs` as `solve` did."""
        T = float(fp.T)
        m = self.cir.n - 1
        rhs = []
        for f in freqs:
            w, _ = pss._forced_replay(fp, f, u_ac, u_points=u_points)
            rhs.append(np.exp(-2j * np.pi * f * T) * np.asarray(w))

        alphas = [np.exp(-2j * np.pi * f * T) for f in freqs]
        tol = max(pss.par.reltol * self.KRYLOV_FACTOR, 1e-14)
        ## ⚠ ON AN OSCILLATOR THE OPERATOR HAS THE ANSWER'S OWN POLE at every
        ## harmonic (see `_check_harmonic`), and a plain solve near one
        ## carries relative error `eta / (2 pi df/f0)`, `eta = |lambda_1 - 1|`
        ## the computed unit multiplier's displacement.  The deflated route
        ## (`_deflated_solve`) borders the pole out and is exact there.
        ## Under the radau default eta ~ 1e-12, so this is correctness
        ## hygiene; the subspace recycling across frequencies is given up on
        ## the autonomous path (one bordered solve per point).
        self.deflated = bool(getattr(pss, 'autonomous', False))
        _dth_f = [None] * len(freqs)
        if self.deflated:
            ## ⚠ ON A STAGED OSCILLATOR the bordered system collapses onto
            ## the total map: with `dtheta = dtheta/dx_0 y_0 + dtheta_f`,
            ## `dtheta_f = -Gt^-1 W f_node` the source's own motion of the
            ## crossings, `(I - a M_tot) y_0 = a (w + P_theta dtheta_f)` --
            ## the deflated solve with the total operator and this source.
            ## Exact against the piecewise-linear forced response (see the
            ## test); the plain deflated solve is 0.3-400x off here.
            _evd = EventColumns.of(pss, fp.width)
            if _evd is not None:
                _Pthd = np.asarray(_evd['P_end'], dtype=complex)
                for i, (f, a) in enumerate(zip(freqs, alphas)):
                    ## (the map's width: gear's pair seeds `(x_0, x_{-1})`)
                    _e0, f_steps = pss._forced_replay(fp, f, u_ac,
                                                      y0=np.zeros(fp.width, dtype=complex),
                                                      collect=True, u_points=u_points)
                    f_nodes = [np.zeros(m, dtype=complex)] + [np.asarray(v_, dtype=complex)[:m]
                                                             for v_ in f_steps]
                    _dth_f[i] = _evd.forced_shift(f_nodes)
                    rhs[i] = np.asarray(rhs[i], dtype=complex) + a * (_Pthd @ _dth_f[i])
            ys = [self._deflated_solve(pss, a, b, transposed=False, tol=tol)
                  for a, b in zip(alphas, rhs)]
            self.matvecs = None
        elif recycle:
            ys, self.matvecs = self._solve_subspace(fp, alphas, rhs, tol)
        else:
            ys, self.matvecs = self._solve_each(fp, alphas, rhs, tol)

        ## one driven replay per frequency turns `y_0` into the period
        ## ⚠ THE BORDERED SIDEBAND RESPONSE.  On a solve whose grid was
        ## landed on state events, a periodic perturbation moves the
        ## crossings: `y_end = M y_0 + w + P_theta dtheta`, and the event
        ## rows close it -- `w_k . y(node_k) = 0` with `y(node) = P_node y_0
        ## + f_node + Pk_node dtheta` (the homogeneous map to the node, the
        ## forced response there, the event column there).  Solved by block
        ## elimination: the m x m solve for the source and for each event
        ## column, then the K x K Schur complement for `dtheta`.  A per-step
        ## saltation reads the dominant multiplier 28 % short of the exact
        ## total (see `_state_event_stage`); the bordered system IS the
        ## linearisation of the solve that produced the orbit.
        _ev = getattr(pss, '_event_columns', None)
        dthetas = [None] * len(freqs)
        if _ev is not None and not self.deflated:
            K = _ev['P_end'].shape[1]
            Wk, nodes = np.asarray(_ev['W']), _ev['nodes']
            for i, (f, a, y0) in enumerate(zip(freqs, alphas, ys)):
                cols = [a * np.asarray(_ev['P_end'][:, k], dtype=complex) for k in range(K)]
                Ycols, _ = self._solve_each(fp, [a] * K, cols, tol)
                _e0, f_steps = pss._forced_replay(fp, f, u_ac, y0=np.zeros(fp.width, dtype=complex),
                                                  collect=True, u_points=u_points)
                f_nodes = [np.zeros(m, dtype=complex)] + [np.asarray(v, dtype=complex)[:m]
                                                         for v in f_steps]
                r = np.zeros(K, dtype=complex)
                S = np.zeros((K, K), dtype=complex)
                for k, nd in enumerate(nodes):
                    Pn = _ev['P_nodes'][nd]
                    _w = Pn.shape[1]          # m on a one-step map, 2m on gear's pair
                    r[k] = Wk[k] @ (Pn @ np.asarray(y0)[:_w] + f_nodes[nd])
                    for l in range(K):
                        S[k, l] = Wk[k] @ (Pn @ np.asarray(Ycols[l])[:_w] + _ev['Pk_nodes'][nd, :, l])
                dth = -np.linalg.solve(S, r)
                dthetas[i] = dth
                ys[i] = np.asarray(y0, dtype=complex) + sum(Ycols[l] * dth[l] for l in range(K))
        elif _ev is not None and self.deflated and _dth_f[0] is not None:
            ## the crossings' motion on the staged oscillator: the state's
            ## part through the total map's sensitivity plus the source's
            _dthx = np.asarray(_ev.dth, dtype=float)
            _w = _dthx.shape[1]           # m on a one-step map, 2m on gear's pair
            for i, y0 in enumerate(ys):
                dthetas[i] = _dthx @ np.asarray(y0, dtype=complex)[:_w] + _dth_f[i]
        ## the crossings' modulation per frequency (fractions of the period
        ## per unit source), None where the solve had no state events
        out = []
        _Pk_fixed = None
        for f, y0, dth in zip(freqs, ys, dthetas):
            _end, ysteps = pss._forced_replay(fp, f, u_ac, y0=y0, collect=True,
                                               u_points=u_points)
            y = np.array([np.asarray(y0)[:m]] + [np.asarray(v)[:m]
                                                 for v in ysteps])
            if dth is not None:
                ## the crossings' motion at every node, AT FIXED TIME --
                ## `Pk_j - xdot_j tau_j^T`: the response of "node j" itself
                ## includes the node's motion along the orbit, O(1) of the
                ## response on a staged oscillator (see
                ## `_fixed_time_event_columns`)
                if _Pk_fixed is None:
                    _Pk_fixed, _t_, _x_ = self._fixed_time_event_columns(pss)
                y = y + np.tensordot(_Pk_fixed[:len(y)], dth, axes=(2, 0))
            out.append(y)
        return out, dthetas

    def adjoint_transfer_row(self, pss, freq, output, recycle_tol=None):
        """Every source to ONE output, in a single transposed solve.

        Returns a row `r` of length `m`: `r[i]` is the small-signal
        response at `output` (at `t = 0`) to a unit source injected at
        reduced coordinate `i` at `freq`. Forward, that is `m` separate
        solves; here it is one.

        ⚠ THIS IS THE ASYMMETRY pnoise IS SHAPED BY, and the reason
        Okumura et al. (1993) reach for the adjoint at all: "it is
        efficient to use the adjoint method ... BECAUSE CIRCUITS HAVE MANY
        NOISE SOURCES." Recycling does not help the forward route, because
        the right-hand side is what changes from source to source.

            output = d^T y_0 = alpha * d^T (I - alpha M)^-1 w(u)
                             = alpha * ((I - alpha M)^-T d)^T W u

        so one transposed solve for `x^a`, then `W^T x^a` from the reverse
        replay, and the whole row falls out.

        The output here is the state at `t = 0`, a single linear
        functional.  A SIDEBAND coefficient `H_l` is a functional
        DISTRIBUTED over the period, whose adjoint takes an injection at
        every step: that is `adjoint_sideband_row`.  Runs under every
        method.

        History: `doc/shooting_history.md`, `PAC.adjoint_transfer_row`.
        """
        import scipy.sparse.linalg as spla
        fp = pss._state_map()
        self._check_circuit(pss)
        self._check_harmonic(pss, freq, 'the adjoint row')
        m = pss.cir.n - 1
        n = fp.width
        alpha = np.exp(-2j * np.pi * float(freq) * float(fp.T))

        d = _output_weights(output, n)

        count = [0]

        def _mv(v):
            count[0] += 1
            return np.asarray(v) - alpha * fp.matvec_transposed(v)

        tol = (self.KRYLOV_FACTOR * pss.par.reltol if recycle_tol is None
               else recycle_tol)
        A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
        ## the same pole as in `solve` and `adjoint_sideband_row`: deflated
        ## on an oscillator, plain (and cheaper) on a driven circuit
        self.deflated = bool(getattr(pss, 'autonomous', False))
        xa = self._pole_solve(pss, A, alpha, d, max(tol, 1e-14),
                              'the adjoint solve')
        self.matvecs = count[0]
        return alpha * pss._forced_replay_transposed(fp, freq, xa)

    def _pole_solve(self, pss, A, alpha, rhs, tol, what):
        """The transposed solve ``(I - alpha M^T) x = rhs`` of an adjoint row:
        DEFLATED on an oscillator, whose operator is singular at every
        harmonic and near-singular around them (the pole is bordered out and
        ``1/(1 - alpha)`` carried analytically), plain GMRES on `A` on a
        driven circuit, where there is no pole and it is both correct and
        cheaper."""
        if getattr(pss, 'autonomous', False):
            return self._deflated_solve(pss, alpha, rhs, transposed=True,
                                        tol=tol)
        return self._gmres_checked(A, rhs, tol, what)

    def adjoint_sideband_row(self, pss, freq, output, sidebands=0):
        """`H_l` rows: every source to ONE output's sideband `l`.

        Returns an array of shape `(len(sidebands), m)`. Entry `[li, i]` is
        the coefficient at sideband `l` of the output at `output`, for a
        unit source injected at reduced coordinate `i` at `freq`:

            H_l = (1/N) sum_n exp(-j l w0 t_n) d^T y_n

        ⚠ THIS IS THE ONE `adjoint_transfer_row` IS NOT.  That row is the
        response at a single instant -- one linear functional, adjointed by
        seeding the reverse pass at the end.  A SIDEBAND is a functional
        DISTRIBUTED over the period, so its adjoint takes an injection at
        EVERY step, and the answer comes in two pieces:

            dH/du  =  [the forced part, from the injected reverse pass]
                    + [alpha * W^T z, with z = (I - alpha M)^-T g]

        where `g` is the reverse pass's own final state -- the sensitivity
        of the functional to the initial state `y_0`, which is itself a
        function of the source through the periodic boundary condition.
        ⚠ Dropping the second term leaves an answer that looks entirely
        reasonable: the two terms are comparable in size (303 against 498
        at `l = 0` on an RC ladder), so neither is a correction to the
        other.

        Still ONE transposed solve per sideband whatever the number of
        sources, which is the property pnoise needs.  Runs under every
        method.

        History: `doc/shooting_history.md`, `PAC.adjoint_sideband_row`.
        """
        import scipy.sparse.linalg as spla
        fp = pss._state_map()

        self._check_circuit(pss)
        self._check_harmonic(pss, freq, 'the sideband row')
        m = pss.cir.n - 1
        n = fp.width
        T = float(fp.T)
        tms = np.asarray(fp.times, dtype=float)
        N = len(fp.steps)
        w0 = 2.0 * np.pi / T
        alpha = np.exp(-2j * np.pi * float(freq) * T)
        ls = np.atleast_1d(np.asarray(sidebands, dtype=int))

        ## ⚠ A HARD BOUND, NOT A HEURISTIC.  Okumura et al. eq. (32): the
        ## maximum frequency the analysis can speak about is the grid's own
        ## `w_max`, so `|l| <= (w_max - w0)/ws`.  You cannot alias down from
        ## above what the grid can represent, and a ratio test on the
        ## accumulated power operates INSIDE this ceiling rather than
        ## instead of it -- an implementation carrying only the ratio test
        ## terminates for the wrong reason.  The grid's ceiling here is its
        ## Nyquist, `N/2` harmonics of the period.
        lmax = N // 2
        bad = ls[np.abs(ls) > lmax]
        if len(bad):
            raise ValueError(
                'PAC: sideband %s is above the grid\'s Nyquist (|l| <= %d '
                'at %d points per period). Nothing can alias down from '
                'above the maximum frequency the grid represents, so this '
                'is not a tolerance to relax -- use a finer period grid.'
                % (bad.tolist(), lmax, N))

        d = _output_weights(output, m)

        count = [0]

        def _mv(v):
            count[0] += 1
            return np.asarray(v) - alpha * fp.matvec_transposed(v)

        A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
        tol = max(self.KRYLOV_FACTOR * pss.par.reltol, 1e-14)

        rows = np.zeros((len(ls), m), dtype=complex)
        for li, l in enumerate(ls):
            ## ⚠ THE PHASE OF THE INPUT COMES OUT FIRST, and getting this
            ## wrong is self-consistent rather than loud.  What is
            ## T-PERIODIC is `v(t) = y(t) exp(-j w t)`, not `y` -- so the
            ## sideband set is the DFT of `v`, which is what `solve` takes.
            ## Decomposing `y` instead gives a Dirichlet kernel smeared
            ## across every `l` whenever `f` is not a multiple of `1/T`,
            ## and it AGREES with a forward reference written the same way,
            ## so only a check against a circuit whose answer is known
            ## independently catches it.
            ## the source couples through every stage/step it reaches (A (x) B
            ## on a coupled tableau), and the output functional is injected at
            ## every node -- one fold for every kind (`_sideband_forced`,
            ## verified vs forward driven solves and the bespoke trbdf2 fold)
            forced, g = pss._sideband_forced(fp, freq, l, d)
            ## ⚠ ON AN OSCILLATOR THIS OPERATOR IS SINGULAR AT EVERY
            ## HARMONIC and near-singular around them, which is exactly
            ## where phase noise is measured.  The deflated route borders
            ## the pole out and carries `1/(1 - alpha)` analytically; on a
            ## driven circuit there is no pole and the plain solve is both
            ## correct and cheaper.
            ## ⚠ ON A STAGED SOLVE THE ROW IS BORDERED: the transpose of
            ## `PAC.solve`'s bordered system, the output read at FIXED times
            ## (`g_theta` over the fixed-time columns), and the event rows'
            ## term -- `-zeta_k W_k` at node k, the source coupling of
            ## `f_node_k` -- as a second reverse pass.  On a driven solve
            ## `z, zeta` come from the block elimination
            ## (`EventColumns.bordered_adjoint`); on an oscillator the system
            ## collapses onto the TOTAL operator, deflated, with `zeta` read
            ## off after it.  Radau/trbdf2 and gear's pair map alike, driven
            ## or free period.  Verified by dual consistency against the
            ## bordered forward solve; unbordered, pnoise on a staged gear
            ## solve is 10-15 % off (driven), and the row missed the forward
            ## solve by 8 % on gear's staged oscillator (2026-09-24, when that
            ## solve first existed).
            _autonomous = getattr(pss, 'autonomous', False)
            _ev = EventColumns.of(pss)
            if _ev is not None and not (fp.is_stage or fp.is_pair
                                        or fp.is_glm):
                _ev = None
            if _ev is not None:
                _wq = pss._period_quadrature(fp)
                cn = np.array([np.exp(-1j * (float(l) * w0 + 2.0 * np.pi * float(freq)) * tms[j])
                               * (1.0 / N if _wq is None else _wq[j]) for j in range(N)])
                _Pkf, _t_, _x_ = self._fixed_time_event_columns(pss)   # the output is read at FIXED times
                g_theta = EventColumns.g_theta(cn, _Pkf, d, N)
                if _autonomous:
                    z = self._deflated_solve(
                        pss, alpha, np.asarray(g, dtype=complex)
                        + np.asarray(_ev.dth, dtype=float).T @ g_theta,
                        transposed=True, tol=tol)
                    zeta = _ev.collapsed_zeta(g_theta, alpha, z)
                else:
                    def _solve_adj(b, k):
                        return self._gmres_checked(
                            A, b, tol, ('the adjoint solve at sideband %d' % l) if k is None
                            else ('the bordered adjoint solve, event %d' % k))
                    z, zeta = _ev.bordered_adjoint(_solve_adj, g, g_theta, alpha)
                forced_ev, _g2 = pss._sideband_forced(
                    fp, freq, l, np.zeros(m), extra=_ev.injection_dict(zeta))
                forced = forced + forced_ev
            else:
                z = self._pole_solve(pss, A, alpha, g, tol,
                                     'the adjoint solve at sideband %d' % l)
            rows[li] = forced + alpha * pss._forced_replay_transposed(
                fp, freq, z)
        self.matvecs = count[0]
        return rows

    ## How small a sideband's contribution must be, relative to the running
    ## total, before the accumulation stops.  Okumura et al.: powers are
    ## "accumulated until their contributions become negligible".
    ALIAS_RATIO_TOL = 1e-9

    def mixer_response(self, pss, f_out, output, sidebands=(-1, 0, 1)):
        """Every input band landing on `f_out` — a `SidebandResponse`.

        For each `l`, the transfer from a source at `f_in = f_out - l f0`
        to the output at `f_out`, which is one
        `adjoint_sideband_row` per sideband because each has its own input
        frequency.  Cost is `len(sidebands)` rows.

        ⚠ A NEGATIVE INPUT BAND IS REFUSED RATHER THAN FOLDED.  For
        `f_out < l f0` the input frequency comes out negative; physically
        that band is the conjugate of `|f_in|`, and quietly taking the
        absolute value would return the right magnitude attached to the
        wrong label -- exactly the mislabelling this object exists to
        prevent.  Ask for the sidebands whose input bands exist.
        """
        self._check_circuit(pss)
        f0 = 1.0 / float(pss.period)
        ls = [int(l) for l in sidebands]
        ins, rows = [], []
        for l in ls:
            f_in = float(f_out) - l * f0
            if f_in < 0.0:
                raise ValueError(
                    'PAC.mixer_response: sideband %d would be fed from '
                    '%.12g Hz, which is negative. That band is the '
                    'conjugate of %.12g Hz, and returning it under the '
                    'label %d would attach a right magnitude to a wrong '
                    'name. Request sidebands whose input bands exist, or '
                    'move f_out.' % (l, f_in, abs(f_in), l))
            row = self.adjoint_sideband_row(pss, f_in, output, sidebands=l)
            ins.append(f_in)
            rows.append(np.asarray(row).reshape(-1))
        return SidebandResponse(f_out, f0, ls, ins, np.array(rows))

    ## ⚠ THE FOLD BELOW IS FOR DRIVEN CIRCUITS.  It is a frequency-conversion
    ## computation and is complete for one; for an AUTONOMOUS oscillator it is
    ## structurally incomplete -- the near-carrier phase-noise skirt is not a
    ## conversion effect (Rizzoli, Mastri & Masotti, MTT 42-807, 1994).  Free-
    ## running phase noise goes through the Floquet/PPV stack instead; see
    ## `oscillator_spectrum` for why the two cannot be unified.
    def pnoise(self, pss, freq, output, ratio_tol=None, maxsidebands=None,
               modulated=False, cyclostationary=False):
        """TIME-AVERAGED output noise PSD at `freq`, sidebands folded in.

            S(f) = sum_l  h_l CY h_l^H ,   h_l = H_l(f - l f0)

        Noise entering at `f - l f0` leaves at `f` through sideband `l`, and
        white sources in disjoint bands are uncorrelated, so the bands add
        in POWER.  Each `h_l` is one adjoint row -- one transposed solve for
        every source in the circuit, which is the whole reason this is
        affordable.

        Returns `(S, sidebands_used)`.  `S` is the one-sided
        **time-averaged** PSD at the output, in the same units as
        `analysis_ss.Noise`'s `Svnout`.  Gated against `analysis_ss.Noise`
        on a linear circuit, where the sidebands vanish and this reduces
        to the stationary answer (Okumura's `p = 1` case).

        The sources' `CY`, three ways:

          * default: STATIONARY SOURCES ONLY, AND IT CHECKS.  A
            bias-dependent `CY` raises (`_cy_reduced`) rather than return a
            number that is quietly the wrong model.
          * `modulated=True`: Hull & Meyer's route -- one stationary source
            at the CYCLE-AVERAGED `CY`, the modulation carried by `H_l`.  It
            keeps the power and drops the correlation between sidebands.
          * `cyclostationary=True`: the construction for a bias-dependent
            `CY`.  A source whose PSD follows the orbit is white noise
            modulated by `B(t) = sqrt(CY(x(t)))`; its band at `f - p f0`
            reaches the output through every modulation harmonic,
            COHERENTLY over `k` and incoherently over `p`, through the SAME
            rows as the stationary fold.  Summing the bands turns the square
            root into the PSD's own harmonics `P_j` (the DFT of `CY(x(t))`,
            no matrix square root, no window count `p`):

                S(f) = sum_{l,l'} a_l P_{l'-l} a_{l'}^H,

            exact on the grid; constant `CY` collapses it to the stationary
            sum.  A coloured source is folded band by band (see
            `_cyclostationary_fold`).  Like the stationary fold it is a
            LOWER bound at a sideband cap.  Coherence is the whole content:
            it differs from `modulated=True` wherever the modulated noise
            crosses a periodically varying transfer, and equals it (only
            `P_0` survives) through a time-invariant one
            (`test_..._cyclostationary_...`).  Cost: white = the cycle
            average's; coloured (any frequency-dependent `CY`, a negligible
            flicker coefficient included) ~6x.
            ⚠ FLICKER: a PSD cannot carry the modulation's SIGN (for a
            coloured source `m xi` and `|m| xi` are different processes), so
            this fold is the `|m|` one: exact when the modulation is
            sign-definite, different physics when it changes sign --
            Okumura's "cannot be modeled as a cyclostationary process by
            using this method, because it has very long time constants".

        ⚠ "TIME-AVERAGED" IS NOT A HEDGE, IT IS THE SPECIFICATION.  Output
        noise is cyclostationary through bias-dependent sources AND through
        the PERIODIC SOURCE-TO-OUTPUT TRANSFER, which applies even when
        every source is stationary.  The time average is sufficient unless
        something downstream tracks the PSD's variation -- a NONLINEAR
        SUBSEQUENT STAGE, or CASCADED STAGES OFF A SHARED REFERENCE
        (Kundert, *Introduction to RF Simulation*).  One number per output
        frequency cannot carry the correlation between frequencies `k f0`
        apart.  (A commercial RF simulator's PNoise computes the same time
        average.)

        ⚠ `maxsidebands` IS AN ACCURACY KNOB HERE AND A REPORTING KNOB IN
        `PAC.solve`.  Noise lives at every frequency, so capping sidebands
        drops power that belonged in the total: `S` becomes a LOWER bound,
        never a cheaper estimate of the same number.

        ⚠ TWO STOPPING RULES, AND THE BOUND IS NOT THE RATIO TEST.  The
        accumulation stops when a sideband pair adds less than `ratio_tol`
        of the running total -- and it can never pass `|l| <= N/2`, the
        grid's own Nyquist (Okumura eq. 32).  `alias_stop` says which
        fired; ending on the bound warns.

        ⚠ HARMONICS: folding puts a copy of a 1/f source's DC singularity
        at every harmonic.  A `freq` on a harmonic where the folded `CY` is
        non-finite or frequency-dependent raises; one just beside it warns.
        Cluster frequencies NEAR each harmonic, never ON it.

        ⚠ AN OSCILLATOR is not this function's problem: its output noise
        is STATIONARY (Demir 2002; `I - M kron M` is singular,
        `test_no_periodic_covariance_exists_for_an_oscillator`), and its
        phase noise is `oscillator_spectrum`'s.  Near a harmonic the
        operator is singular with the PPV as its null vector, which a plain
        solve shows as "flat PSD curves or curves with unexpected slope
        near the oscillation frequency" (Gourary et al.).  The rows use the
        deflated solve (`_deflated_solve`, Gourary et al. eq. 27/28,
        bordered with BOTH null vectors); `PAC.deflated` says which route
        ran.

        History: `doc/shooting_history.md`, `PAC.pnoise`.
        """
        self._check_circuit(pss)
        ## pnoise folds sidebands through the ADJOINT (adjoint_sideband_row ->
        ## _forced_replay_transposed), whose two-stage chained transpose is
        ## not built for TR-BDF2, so it falls back to a Gear-2 twin -- see
        ## `_adjoint_host`.  (covariance/oscillator_covariance use the built
        ## TR-BDF2 Lyapunov injection via `_lyapunov_host`.)
        pss = pss._adjoint_host()
        fp = pss.factored_period()
        m = pss.cir.n - 1
        N = len(fp.steps)
        T = float(fp.T)
        f0 = 1.0 / T
        tol = self.ALIAS_RATIO_TOL if ratio_tol is None else float(ratio_tol)
        lmax = N // 2 if maxsidebands is None else min(int(maxsidebands),
                                                       N // 2)

        w = 2.0 * np.pi * float(freq)
        ## ⚠ `modulated=True` IS HULL & MEYER'S ROUTE, NOT A TOLERANCE
        ## RELAXATION.  Off, a bias-dependent `CY` raises, because the
        ## stationary sum would be the wrong model.  On, the source is
        ## replaced by ONE stationary source at the CYCLE-AVERAGED bias and
        ## the modulation is carried by `H_l` -- which is the standard
        ## treatment of exactly this case, and the only route to MOS
        ## pnoise, since no physically correct MOS noise model has a
        ## state-independent `CY`.
        colour = None
        if cyclostationary:
            ## the stop rule and the harmonic probes below run on the
            ## cycle-averaged power (the modulation's B_0 B_0^H); the fold
            ## itself is the convolution after the rows are gathered.  The
            ## colour model is fitted ONCE here and serves both, per element
            ## so independent sources ADD in the coloured fold (see
            ## `_cy_components_model`); its call is the summed model the
            ## stop rule reads
            colour = self._cy_components_model(pss, float(freq), f0)
            if colour is None:
                cyfn = self._cy_cycle_averaged
            else:
                fp_ = pss.factored_period()
                hs_ = np.diff(np.asarray(fp_.times, dtype=float))
                def cyfn(pss_, w_, _m=colour, _h=hs_):
                    Cs = _m(w_)
                    ns = min(len(_h), Cs.shape[0])
                    return np.einsum('k,kij->ij', _h[:ns], Cs[:ns]) / float(_h[:ns].sum())
        else:
            cyfn = (self._cy_cycle_averaged if modulated else self._cy_reduced)
        cy = cyfn(pss, w)

        ## ⚠⚠ ON A HARMONIC, A SIDEBAND FOLDS THE SOURCES TO DC -- AND
        ## SOME DEVICE MODELS ARE NOT DEFINED THERE.  Sideband `l`
        ## evaluates `CY` at `f - l f0`, so `f = k f0` evaluates it at
        ## ZERO: a `1/f` term is infinite there, and a flicker term with its
        ## coefficient set to ZERO is `0/0` = `nan` (`PspMosLongChannel`,
        ## `fnt=1, nfa=0`) -- a caller who sets `nfa = 0` believing flicker
        ## is off still gets `nan`.  White sources fold to DC harmlessly, so
        ## this checks the SOURCES at the frequency actually used rather
        ## than refusing a harmonic on principle.
        f0_ = 1.0 / float(pss.period)
        lscan = max(1, int(maxsidebands or 8))
        offs = np.abs(float(freq) - np.arange(-lscan, lscan + 1) * f0_)
        near = float(np.min(offs))
        if near <= self.HARMONIC_GUARD * f0_:
            probe = cyfn(pss, 2.0 * np.pi * near)
            if not np.all(np.isfinite(np.asarray(probe))):
                raise ValueError(
                    'PAC.pnoise: %.12g Hz sits on a harmonic of %.12g Hz, '
                    'so a sideband folds the noise sources to DC -- and at '
                    'DC this circuit\'s CY is not finite. A 1/f term is '
                    'infinite there; a flicker term whose COEFFICIENT IS '
                    'ZERO is 0/0 and gives nan, so disabling flicker does '
                    'not avoid this. Offset from the harmonic: a commercial RF simulator\'s '
                    'own advice is to cluster frequencies NEAR each '
                    'harmonic and never place one ON it.'
                    % (float(freq), f0_))
            ## ⚠⚠ A FINITE PROBE IS NOT A SAFE ONE.  `1/T` rounds, so at
            ## `f = f0` the folded band sits ~1e-11 Hz from DC, not ON it: a
            ## 1/f source there is finite and enormous (6.3e-2 V^2/Hz against
            ## 9.2e-15 at 0.1 % either side).  So a frequency-DEPENDENT CY at
            ## the folded band refuses too; a white one stays allowed.
            probe_ref = cyfn(pss, 2.0 * np.pi * max(abs(float(freq)) * 2.0, f0_))
            if not np.allclose(np.asarray(probe), np.asarray(probe_ref),
                               rtol=1e-9, atol=0.0):
                raise ValueError(
                    'PAC.pnoise: %.12g Hz sits on harmonic %d of %.12g Hz, so '
                    'a sideband folds the noise sources to %.3g Hz -- DC up to '
                    'rounding -- and this circuit\'s CY is frequency-dependent '
                    'there: a 1/f source is read next to its singularity and '
                    'the fold returns a finite, absurd number (measured '
                    '6.3e-2 V^2/Hz against 9.2e-15 at 0.1 %% either side). '
                    'Offset from the harmonic, or use PAC.sampled_variance, '
                    'whose explicit fmin keeps every band off DC.'
                    % (float(freq), int(round(abs(float(freq)) / f0_)), f0_,
                       near))

        ## ⚠ THE STEEP REGION BESIDE A HARMONIC IS A SWEEP HAZARD RATHER
        ## THAN A WRONG NUMBER, so it warns instead of raising: the VALUE is
        ## right (2 % above the plateau at `f0 + 0.01` Hz with a real flicker
        ## source), but a grid that lands there by accident integrates a
        ## spike it never resolved.
        elif near < 1e-6 * f0_ and float(freq) > 0.0:
            cy_hi = cyfn(pss, 2.0 * np.pi * max(float(freq) * 2.0, f0_))
            if not np.allclose(cy, cy_hi, rtol=1e-9, atol=0.0):
                warnings.warn(
                    'PAC.pnoise: %.12g Hz is %.3g Hz from a harmonic of '
                    '%.12g Hz and a source has a frequency-dependent CY, '
                    'so the folded density varies steeply here. The VALUE '
                    'is correct; a swept grid landing this close will '
                    'misrepresent the integrated total. Cluster near each '
                    'harmonic deliberately rather than by accident.'
                    % (float(freq), near, f0_), RuntimeWarning, stacklevel=2)

        total = 0.0
        used = []
        quiet = 0
        self.alias_stop = 'bound'
        rows = {}
        for l in range(0, lmax + 1):
            step = 0.0
            for sl in ((0,) if l == 0 else (l, -l)):
                fin = float(freq) - sl * f0
                h = self.adjoint_sideband_row(pss, fin, output, sl)[0]
                rows[sl] = np.asarray(h, dtype=complex)
                step += float(np.real(h @ cyfn(
                    pss, 2.0 * np.pi * fin) @ np.conj(h)))
                used.append(sl)
            total += step
            if total > 0 and abs(step) < tol * abs(total):
                ## ⚠ TWO QUIET PAIRS, NOT ONE.  A single sideband can come
                ## back near zero by symmetry while its neighbours do not,
                ## and stopping there would truncate a series that had not
                ## converged.
                quiet += 1
                if quiet >= 2:
                    self.alias_stop = 'ratio'
                    break
            else:
                quiet = 0
        self.sidebands_used = used
        if cyclostationary:
            total = self._cyclostationary_fold(pss, float(freq), rows, model=colour)
        ## ⚠ WHICH RULE STOPPED IT IS PART OF THE ANSWER.  Ending on the
        ## ratio test means the series converged; ending on the Nyquist
        ## bound means the grid ran out before the series did, and the
        ## number is a LOWER bound on the folded noise -- every sideband
        ## above the grid's own maximum frequency is missing, not small.
        ## A strongly switching circuit does this readily.
        if self.alias_stop == 'bound' and lmax > 0:
            warnings.warn(
                'PAC.pnoise: the sideband accumulation stopped at the '
                "grid's Nyquist (|l| = %d at %d points per period), not "
                'because the contributions became negligible. Sidebands '
                'above the grid\'s maximum frequency are MISSING rather '
                'than small, so this is a lower bound on the folded noise. '
                'Re-solve the PSS on a finer period grid and compare.'
                % (lmax, N),
                RuntimeWarning, stacklevel=2)
        return total, used

    def _cy_harmonics(self, pss, w):
        """`P_j`: the Fourier coefficient matrices of `CY(x(t), w)` over the
        orbit, `(N, n, n)` indexed like `numpy.fft.fftfreq`.  `P_0` is the
        cycle average; `P_j` with `j != 0` carry the modulation and vanish
        for a bias-independent source.  No square root: the fold uses the
        PSD's own harmonics (`a P a^H`), exact on the grid at any window.
        (A sqrt-modulation route is not: the root of a PSD that crosses
        zero has a kink and a slowly decaying harmonic tail, which the
        sideband window truncates.)

        History: `doc/shooting_history.md`, `PAC._cy_harmonics`."""
        return self._period_dft(pss, self._cy_samples(pss, w))

    def _cy_samples(self, pss, w):
        """`CY(x(t_k), w)` over the orbit, reduced, `(N, n, n)` complex."""
        fp = pss.factored_period()
        irn = pss.irefnode
        xs = np.asarray(pss.waveform[1], dtype=float)
        nsamp = len(fp.steps)
        Cs = []
        for k in range(nsamp):
            xr = np.asarray(xs[:, k], dtype=float).ravel()
            xf = xr if xr.shape[0] == pss.cir.n else np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            cyk = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cyk,) = remove_row_col((cyk,), irn, pss.toolkit)
            Cs.append(np.asarray(cyk, dtype=complex))
        return np.asarray(Cs, dtype=complex)

    def _period_dft(self, pss, S):
        """Fourier coefficients over the period of samples `S` `(N, ...)` taken
        at `fp.times[:N]`, laid out like `numpy.fft.fftfreq`.  Uniform grid:
        the index DFT, unchanged.  Non-uniform: the weighted sum at the TRUE
        times with `PSS._period_quadrature`'s weights -- O(N^2), paid only by
        a caller who chose that grid."""
        S = np.asarray(S, dtype=complex)
        fp = pss.factored_period()
        wq = pss._period_quadrature(fp)
        if wq is None or S.shape[0] != len(wq):
            return np.fft.fft(S, axis=0) / S.shape[0]
        N = S.shape[0]
        tms = np.asarray(fp.times, dtype=float)
        ks = np.fft.fftfreq(N, d=1.0 / N)
        E = np.exp(-2j * np.pi * np.outer(ks, (tms[:N] - tms[0]) / float(tms[N] - tms[0]))) * wq[None, :]
        return np.tensordot(E, S, axes=(1, 0))

    @staticmethod
    def _sqrt_harmonics_of(Cs, dft=None):
        """The DFT of the symmetric square root of the sampled `CY` (see
        `_cy_sqrt_harmonics`), for a `(N, n, n)` array already in hand."""
        Bs = []
        for cyk in Cs:
            cyk = 0.5 * (cyk + cyk.conj().T)
            lam, U = np.linalg.eigh(cyk)
            lam = np.clip(np.real(lam), 0.0, None)
            Bs.append((U * np.sqrt(lam)[None, :]) @ U.conj().T)
        Bs = np.asarray(Bs, dtype=complex)
        return (np.fft.fft(Bs, axis=0) / Bs.shape[0]) if dft is None else dft(Bs)

    def _cy_colour_model(self, pss, f, f0):
        """Fit `CY(x(t), w) = A(t) + B(t) (w1/w)^ef` entry by entry from three
        frequencies and verify at a fourth; return a callable `w -> (N, n, n)`
        or None when the fit fails anywhere (the caller then evaluates the
        circuit per band, as before).  The exponent is per entry, found by
        a bracketed root find on the ratio of differences, so a mix of
        flicker exponents across sources is fine; a white entry (B = 0)
        needs no exponent."""
        ws = self._colour_fit_frequencies(f, f0)
        fit = self._colour_fit([self._cy_samples(pss, w) for w in ws], ws)
        if fit is None:
            return None
        A, B, EF = fit
        w1 = ws[0]
        def model(w):
            ## |w| -- see `_cy_components_model`: a negative band frequency
            ## otherwise gives NaN and silently disables pnoise's ratio stop
            with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
                return A + B * (np.float64(w1) / np.abs(np.float64(w))) ** EF
        return model

    @staticmethod
    def _colour_fit_frequencies(f, f0):
        ## three to fit, two to verify: one BETWEEN the fit points and one
        ## at the FAR end of the band range the fold reaches (up to
        ## ~(N/2 + lmax) f0), so a shape that is not thermal-plus-flicker
        ## is caught where the model would have been extrapolating
        f = abs(float(f))
        return [2.0 * np.pi * x for x in (max(f, 1e-3 * f0), 3.0 * f0 + f,
                                          10.0 * f0 + f, 2.0 * f0 + f,
                                          150.0 * f0 + f)]

    @staticmethod
    def _colour_fit(Cs, ws):
        """`(A, B, EF)` with `C(w) = A + B (w1/w)^EF` entry by entry, from
        `Cs` = samples `(N, n, n)` at the five `ws` of
        `_colour_fit_frequencies` (three fit, two verify); None when the
        shape is not thermal-plus-power-law anywhere."""
        from scipy.optimize import brentq
        C1, C2, C3, C4, C5 = Cs
        N, n, _ = C1.shape
        A = np.zeros_like(C1); B = np.zeros_like(C1); EF = np.zeros((N, n, n))
        scale = max(float(np.max(np.abs(C1))), 1e-300)
        w1, w2, w3, w4, w5 = ws
        for k in range(N):
            for i in range(n):
                for j in range(n):
                    c1, c2, c3 = C1[k, i, j], C2[k, i, j], C3[k, i, j]
                    if abs(c1 - c2) <= 1e-12 * scale and abs(c2 - c3) <= 1e-12 * scale:
                        A[k, i, j] = c1
                        continue
                    ratio = (c1 - c2) / (c2 - c3)
                    def g(ef, ratio=ratio):
                        g1, g2, g3 = 1.0, (w1 / w2) ** ef, (w1 / w3) ** ef
                        return float(np.real((g1 - g2) / (g2 - g3) - ratio))
                    try:
                        ef = brentq(g, 0.05, 4.0, xtol=1e-12)
                    except ValueError:
                        return None
                    g2 = (w1 / w2) ** ef
                    Bv = (c1 - c2) / (1.0 - g2)
                    A[k, i, j] = c1 - Bv; B[k, i, j] = Bv; EF[k, i, j] = ef
        for wv, Cv in ((w4, C4), (w5, C5)):
            if float(np.max(np.abs(A + B * (w1 / wv) ** EF - Cv))) > 1e-8 * scale:
                return None
        return A, B, EF

    @classmethod
    def _leaf_cy_stamps(cls, cir, x, w, prefix=()):
        """Yield `(key, G)`: each LEAF element's `CY(x, w)` stamped into
        `cir`'s full `n x n` space, recursing into sub-circuits.  Their sum
        is `cir.CY(x, w)` (the elements are independent by that method's
        own contract)."""
        n = cir.n
        idx = cir._map_indices_2d
        for inst, el in cir.elements.items():
            rc = idx.get(inst)
            if rc is None:
                continue
            rows, cols = rc
            subx = np.asarray(x)[cir.elementnodemap[inst]]
            if getattr(el, 'elements', None):
                for key, Gc in cls._leaf_cy_stamps(el, subx, w, prefix + (inst,)):
                    G = np.zeros((n, n), dtype=complex)
                    np.add.at(G, (rows, cols), np.asarray(Gc).ravel())
                    yield key, G
            else:
                G = np.zeros((n, n), dtype=complex)
                np.add.at(G, (rows, cols),
                          np.asarray(el.CY(subx, w), dtype=complex).ravel())
                yield prefix + (inst,), G

    @classmethod
    def _leaf_noise_amplitudes(cls, cir, x, w, prefix=()):
        """Yield `(key, W)`: each leaf element's SIGNED coloured-noise
        amplitudes (`Element.noise_amplitudes`, where it has them) in `cir`'s
        full `n`-row space, `(n, S)`; keyed like `_leaf_cy_stamps`."""
        n = cir.n
        for inst, el in cir.elements.items():
            if cir._map_indices_2d.get(inst) is None:
                continue
            nodemap = np.asarray(cir.elementnodemap[inst])
            subx = np.asarray(x)[nodemap]
            if getattr(el, 'elements', None):
                inner = cls._leaf_noise_amplitudes(el, subx, w, prefix + (inst,))
            else:
                fn = getattr(el, 'noise_amplitudes', None)
                Wc = fn(subx, w) if fn is not None else None
                inner = [] if Wc is None else [(prefix + (inst,), Wc)]
            for key, Wc in inner:
                Wc = np.asarray(Wc, dtype=complex)
                W = np.zeros((n, Wc.shape[1]), dtype=complex)
                np.add.at(W, nodemap, Wc)
                yield key, W

    def _signed_amplitudes(self, pss, w1, flicker, states=None):
        """`{key: (K, m, S)}` for the flicker components whose element states
        its signed amplitudes AND whose amplitudes rebuild the component:
        `W W^dagger = B` at the fit frequency, to 1e-6.  Anything else keeps
        the square root of its PSD -- the |m| process, warned on as before."""
        irn = pss.irefnode
        keep = np.array([i for i in range(pss.cir.n) if i != irn])
        acc = {}
        for xf in self._orbit_states(pss, states):
            for key, W in self._leaf_noise_amplitudes(pss.cir, xf, w1):
                acc.setdefault(key, []).append(W[keep])
        out = {}
        for key, B, _EF in flicker:
            if key not in acc:
                continue
            W = np.asarray(acc[key], dtype=complex)
            if W.shape[0] != B.shape[0] or not np.all(np.isfinite(W)):
                continue
            rebuilt = np.einsum('kis,kjs->kij', W, W.conj())
            if float(np.max(np.abs(rebuilt - B))) <= 1e-6 * float(np.max(np.abs(B))):
                out[key] = W
        return out

    @staticmethod
    def _warn_signed_unused(model, where):
        """Warn when an element STATED its signed amplitudes and the fold
        factors that component by sqrt(PSD) anyway: a silent fallback
        reproduces the sign-blind answer, which looks like agreement.

        History: `doc/shooting_history.md`, `PAC._warn_signed_unused`."""
        signed = getattr(model, 'amplitude', None) or {}
        lost = [key for key, B, EF in (getattr(model, 'flicker', None) or [])
                if key in signed and PAC._uniform_exponent(B, EF) is None]
        if lost:
            warnings.warn(
                '%s: %s states SIGNED coloured-noise amplitudes, but its '
                'power-law exponent is not uniform across its entries, so the '
                'component is evaluated per band from sqrt(PSD) -- the SIGN-'
                'BLIND fold (the |m| process).  The result is the pre-2026-09-19 '
                'one for this component, not the signed physics.'
                % (where, ', '.join('.'.join(k) for k in lost)),
                RuntimeWarning, stacklevel=4)

    @staticmethod
    def _orbit_states(pss, states=None):
        """Full-width state vectors to sample `CY` at: the stored orbit
        `x(t_k)`, `k = 0..N-1` (the default), or the given `states`."""
        n = pss.cir.n
        irn = pss.irefnode
        if states is None:
            xs = np.asarray(pss.waveform[1], dtype=float)
            nsamp = len(pss.factored_period().steps)
            states = [xs[:, k] for k in range(nsamp)]
        out = []
        for xr in states:
            xr = np.asarray(xr, dtype=float).ravel()
            out.append(xr if xr.shape[0] == n
                       else np.concatenate((xr[:irn], np.zeros(1), xr[irn:])))
        return out

    def _element_cy_samples(self, pss, w, states=None):
        """`{key: (K, m, m)}` -- each leaf element's reduced `CY(x, w)` at the
        orbit samples `x(t_k)` (the default, indexed like `_cy_samples`, which
        they sum to) or at the given `states`."""
        irn = pss.irefnode
        n = pss.cir.n
        keep = np.array([i for i in range(n) if i != irn])
        out = {}
        for xf in self._orbit_states(pss, states):
            for key, G in self._leaf_cy_stamps(pss.cir, xf, w):
                out.setdefault(key, []).append(G[np.ix_(keep, keep)])
        return {key: np.asarray(v, dtype=complex) for key, v in out.items()}

    def _cy_at_states(self, pss, w, states=None):
        """The whole circuit's reduced `CY(x, w)` at the orbit samples or at
        `states`, `(K, m, m)` -- `_cy_samples` generalised."""
        irn = pss.irefnode
        Cs = []
        for xf in self._orbit_states(pss, states):
            cyk = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cyk,) = remove_row_col((cyk,), irn, pss.toolkit)
            Cs.append(np.asarray(cyk, dtype=complex))
        return np.asarray(Cs, dtype=complex)

    def _cy_components_model(self, pss, f, f0, states=None):
        """`CY(x(t), w)` as INDEPENDENT components -- one white and one
        coloured part per leaf element -- for the coloured folds.

        ⚠⚠ Independent sources whose modulations differ do not add under
        a joint square root: `sqrt(A(t) + B)` cross-couples them at
        `(t, t')` (a switch's white noise plus a 1/f source at one node:
        +7.3 % of the total).  One root per component restores
        additivity.

        Returns a callable `w -> (N, m, m)` (the summed `CY`, the contract
        `_cy_colour_model` had) carrying `white` (summed), `white_parts`
        `[(key, A)]`, `flicker` `[(key, B, EF)]` (`C = A + B (w1/w)^EF` per
        element, fitted and verified as `_colour_fit`), `perband` (keys
        whose colour did not fit: evaluated per band) and `w1`.  ⚠ Within
        one element all white terms share a root, as do all power-law terms
        -- independence is resolved to element x {white, coloured}.

        History: `doc/shooting_history.md`, `PAC._cy_components_model`.
        """
        ws = self._colour_fit_frequencies(f, f0)
        per_w = [self._element_cy_samples(pss, w, states) for w in ws]
        keys = [key for key in per_w[0]
                if any(np.any(per_w[i][key]) for i in range(len(ws)))]
        m = pss.cir.n - 1
        N = (len(pss.factored_period().steps) if states is None
             else len(states))
        ## ⚠ THE ELEMENTS MUST BE THE CIRCUIT'S CY.  A circuit whose `CY` is
        ## not the sum of its leaf elements' (an override, a batched toolkit
        ## group) cannot be split, and splitting it anyway would silently
        ## analyse a different noise model -- so check at two fit frequencies
        ## and fall back to the whole-circuit model, saying why.
        for i in (0, len(ws) - 1):
            whole = self._cy_at_states(pss, ws[i], states)
            parts = sum((per_w[i][key] for key in per_w[i]),
                        np.zeros_like(whole))
            scale = max(float(np.max(np.abs(whole))), 1e-300)
            if float(np.max(np.abs(parts - whole))) > 1e-9 * scale:
                warnings.warn(
                    'PAC: this circuit\'s CY is not the sum of its elements\' '
                    '(an override or a batched group?), so the coloured fold '
                    'cannot take one square root per independent source and '
                    'uses ONE root of the summed CY -- independent sources '
                    'with different modulations are then not additive '
                    '(measured +7.3 % of the total on a switch + 1/f source).',
                    RuntimeWarning, stacklevel=3)
                if states is not None:
                    return None
                return self._cy_colour_model(pss, f, f0)
        white = np.zeros((N, m, m), dtype=complex)
        white_parts, flicker, perband = [], [], []
        for key in keys:
            fit = self._colour_fit([per_w[i][key] for i in range(len(ws))], ws)
            if fit is None:
                perband.append(key)
                continue
            A, B, EF = fit
            if np.any(A):
                white = white + A
                white_parts.append((key, A))
            if np.any(B):
                flicker.append((key, B, EF))
        if perband:
            warnings.warn(
                'PAC: the noise of %s is not thermal-plus-power-law, so it is '
                'evaluated per band with ONE square root per element: '
                'independent sources INSIDE such an element are not split '
                '(measured 4.2e-4 on an EKV stage, thermal + flicker).'
                % ', '.join('.'.join(k) for k in perband),
                RuntimeWarning, stacklevel=3)
        w1 = ws[0]

        def model(w):
            tot = white.copy()
            ## ⚠ |w|: the stop rule asks at NEGATIVE band frequencies
            ## (f - l f0 < 0), where a non-integer power of a negative base is
            ## NaN and the ratio stop never fires
            ## ⚠ and numpy division: exactly ON a harmonic the folded band is
            ## w = 0, where a Python float division raises ZeroDivisionError
            ## from inside pnoise's harmonic guard instead of letting it see
            ## the non-finite CY and refuse by name
            with np.errstate(divide='ignore', invalid='ignore', over='ignore'):
                ratio = np.float64(w1) / np.abs(np.float64(w))
                for _key, B, EF in flicker:
                    tot = tot + B * ratio ** EF
            if perband:
                ew = self._element_cy_samples(pss, w, states)
                for key in perband:
                    tot = tot + ew[key]
            return tot
        model.white = white
        model.white_parts = white_parts
        model.flicker = flicker
        model.perband = perband
        ## ⚠ THE SIGN: where the element states its coloured AMPLITUDES, the
        ## folds factor the component with them instead of with sqrt(PSD) --
        ## see `Element.noise_amplitudes` (hdl.py)
        model.amplitude = self._signed_amplitudes(pss, w1, flicker, states)
        model.w1 = w1
        return model

    @staticmethod
    def _uniform_exponent(B, EF):
        """The one power-law exponent of a component, or None when its
        non-zero entries carry different ones (then `sqrt(B (w1/w)^EF)` is
        not `(w1/w)^(ef/2) sqrt(B)` and must be taken per band)."""
        ## ⚠ ONLY ENTRIES THAT CARRY WEIGHT VOTE, weighted by what a wrong
        ## exponent COSTS.  The exponent is fitted from differences of `CY`,
        ## so a tiny entry (a MOS flicker source at the sample where Vds
        ## crosses zero) has its exponent in the rounding of the white part
        ## beside it, and that noise goes as 1/weight: no weight cut-off
        ## separates them.  Giving entry i the exponent `ref` misstates the
        ## component by `r_i |(w1/w)^d_i - 1| ~ r_i d_i |ln(w1/w)|` of its
        ## scale (`r` relative weight, `d` deviation); bounded over 50
        ## e-folds of band frequency and held to 1e-9, a genuinely different
        ## exponent (d ~ 1) still fails from a weight of 2e-11 up.  The
        ## reference is the LARGEST entry's exponent, not the first's.
        ## History: `doc/shooting_history.md`, `PAC._uniform_exponent`.
        aB = np.abs(B)
        if not np.any(aB > 0):
            return 0.0
        ref = float(np.real(EF.flat[int(np.argmax(aB))]))
        cost = 50.0 * (aB / float(aB.max())) * np.abs(EF - ref)
        return ref if float(np.max(cost)) <= 1e-9 else None

    def _cy_sqrt_harmonics(self, pss, w):
        """`B_k`: the DFT of the symmetric square root of `CY(x(t), w)` over
        the orbit, `(N, n, n)`, for the band-resolved (coloured) fold."""
        Bs = []
        for cyk in self._cy_samples(pss, w):
            cyk = 0.5 * (cyk + cyk.conj().T)
            lam, U = np.linalg.eigh(cyk)
            lam = np.clip(np.real(lam), 0.0, None)
            Bs.append((U * np.sqrt(lam)[None, :]) @ U.conj().T)
        Bs = np.asarray(Bs, dtype=complex)
        return self._period_dft(pss, Bs)

    def _cyclostationary_fold(self, pss, freq, rows, model=None):
        """`S(f) = sum_{l,l'} a_l Q_{l,l'} a_{l'}^H` over the gathered
        sideband rows (`rows[l]` = the row for a source at `f - l f0`,
        output at `f`).  WHITE source: `Q_{l,l'} = P_{l'-l}`, the DFT of
        `CY(x(t))` itself -- no square root, exact on the grid.  COLOURED
        source (`CY` depends on `w`; detected by comparing two bands): the
        white band `p = l + k` shared by rows `l` and `l'` carries its OWN
        `CY`, so `Q_{l,l'} = sum_k B_k^{(l+k)} B_{k+l-l'}^{(l+k) H}` with
        `B^{(p)}` the sqrt-DFT at the band's frequency `|f - p f0|`, summed
        over ALL `N` modulation harmonics `k` (which is what makes the
        square root exact here; a window on `k` is not).  ⚠ On a flicker
        source `||P_0||` differs 24x across the bands the fold sums, so
        "the band of l" is not an approximation to use; the band-resolved
        form is pinned against the stationary fold of a stationary FLICKER
        source through the same multiplier.

        The colour model (`_cy_components_model`, fitted once in `pnoise`
        and shared with the stop rule) and a vectorised pair sum keep the
        coloured call near the white one, exact to 1e-11 against the
        per-band evaluation, which remains the fallback for a colour the
        model does not fit.

        History: `doc/shooting_history.md`, `PAC._cyclostationary_fold`."""
        f0 = 1.0 / float(pss.period)
        ls = sorted(rows)
        f = float(freq)
        ## coloured or white?  two bands, same test the stationary path uses
        w_a = 2.0 * np.pi * abs(f - ls[0] * f0)
        w_b = 2.0 * np.pi * max(abs(f) * 2.0, f0)
        Pa = self._cy_harmonics(pss, w_a)
        Pb = self._cy_harmonics(pss, w_b)
        coloured = not np.allclose(Pa, Pb, rtol=1e-9, atol=0.0)
        total = 0.0
        if not coloured:
            P = Pa
            N = P.shape[0]
            for l in ls:
                for lp in ls:
                    total += complex(rows[l] @ P[(lp - l) % N] @ np.conj(rows[lp]))
            return float(np.real(total))
        ## band-resolved: every white band the modulation harmonics reach.
        ## The cost is the circuit's CY, not the algebra: every colour in
        ## the library is thermal plus flicker in 1/f^ef, so the model fixes
        ## each entry's shape from three evaluations per sample (A + B
        ## (w1/w)^ef, ef by a root find), verifies it at further
        ## frequencies, and gives every band with no further circuit calls;
        ## a colour not of that shape gets the full evaluation.
        Nn = Pa.shape[0]
        if model is None:
            model = self._cy_components_model(pss, f, f0)
        ## ⚠ A SPECIFICATION LIMIT, NOT AN IMPLEMENTATION ONE: a coloured
        ## source under a modulation that CHANGES SIGN is not representable
        ## by any fold built from a PSD -- R(t,t') = m(t) m(t') R_c(t-t')
        ## keeps the sign product and CY cannot carry it -- so for a source
        ## that states no signed amplitudes this fold computes the |m|
        ## process.  A device's own 1/f current is a conductance fluctuation
        ## TIMES the current and follows its sign.  Where the element states
        ## its signed amplitudes (`Element.noise_amplitudes`) the folds use
        ## them and nothing below applies.  (White sources are untouched:
        ## uncorrelated across the period, no sign product survives.)
        ## The sign is invisible here; its NECESSARY condition is a PSD that
        ## touches zero along the orbit with a KINK in its square root, so
        ## that is warned on.  The touch threshold: a zero crossing SAMPLED
        ## on an N-point grid bottoms out near (pi/N)^2 of the maximum, a
        ## sign-definite PSD with a ten-fold swing sits at 1e-2 -- so 1e-2;
        ## a heuristic, and a warning for that reason.
        ## ⚠ NOT WHEN EVERY COLOURED COMPONENT CARRIES ITS SIGN: then nothing
        ## below takes a square root of a PSD and there is nothing to warn of
        _signed = getattr(model, 'amplitude', None) or {}
        self._warn_signed_unused(model, 'PAC.pnoise(cyclostationary=True)')
        _all_signed = (getattr(model, 'flicker', None) is not None
                       and not model.perband
                       and all(k_ in _signed and self._uniform_exponent(B_, E_) is not None
                               for k_, B_, E_ in model.flicker))
        Cs0 = np.asarray([np.abs(np.diag(np.fft.ifft(Pa, axis=0)[k])) for k in range(Nn)])
        dmax = Cs0.max(axis=0)
        touches = (dmax > 0) & (Cs0.min(axis=0) <= 1e-2 * dmax)
        ## ⚠ THE ORDER OF THE ZERO: a LINEAR sign crossing gives sqrt(PSD) a
        ## first-derivative KINK, a sign-definite quadratic touch a smooth
        ## minimum.  The circular second difference of sqrt(PSD) divided by
        ## h/T and by the maximum is a DERIVATIVE JUMP: grid-independent at a
        ## kink (~4 pi for a sinusoidal slope) and falling as h/T where
        ## smooth, so 3 separates them down to ~50 points per period (a raw
        ## threshold would encode the grid).  ⚠ STILL NECESSARY, NOT
        ## SUFFICIENT, AND THE SENSITIVITY RUNS INVERSE TO THE EFFECT: a
        ## crossing flatter than linear (an LO shaped v |v|^(p-1), p > 1)
        ## stays O(1) wrong while the indicator falls by orders.  A quiet
        ## warning is not evidence of a small discrepancy.
        kinked = np.zeros_like(touches)
        hT = 1.0 / float(Nn)
        for jj in np.where(touches)[0]:
            sq = np.sqrt(Cs0[:, jj])
            d2 = np.abs(sq - 0.5 * (np.roll(sq, 1) + np.roll(sq, -1)))
            kinked[jj] = bool(d2.max() / (sq.max() * hT) > 3.0)
        if bool(np.any(kinked)) and not _all_signed:
            warnings.warn(
                'PAC.pnoise(cyclostationary=True): a COLOURED source whose PSD '
                'touches zero along the orbit -- if its modulation changes sign '
                '(a switching gain), no PSD-specified model can represent the '
                'coloured process (Okumura eq. 23 in concrete form), and this '
                'fold computes the |m| one (its square root has a first-derivative '
                'kink at the zero, the signature of a LINEAR sign crossing; a '
                'necessary condition -- a shallow crossing shows no kink and '
                'errs MORE): measured 0.56x and 1.33x of the signed '
                'physics at two offsets on a flicker source through a '
                'zero-crossing gain -- EITHER direction, the sign of the '
                'discrepancy is set by the offset, not the mechanism -- and '
                'exact for a sign-definite one. Only the element knows the sign.',
                RuntimeWarning, stacklevel=3)
        ks = np.fft.fftfreq(Nn, d=1.0 / Nn).astype(int)
        pmin = min(ls) + int(ks.min()); pmax = max(ls) + int(ks.max())
        wband = lambda p: 2.0 * np.pi * abs(f - p * f0)
        ## ⚠⚠ ONE SQUARE ROOT PER INDEPENDENT COMPONENT, NOT OF THE SUM: a
        ## joint `sqrt(CY)` makes independent sources with different
        ## modulations NON-ADDITIVE -- see `_cy_components_model`.
        ## The white parts stay in the exact P-form (linear in `CY`, so
        ## additive already); each coloured part gets its own root, scaled
        ## per band when its exponent is uniform (`sqrt(c B) = sqrt(c)
        ## sqrt(B)`, no per-band eigendecomposition).
        ## the period's Fourier coefficients: the index DFT on a uniform grid,
        ## the trapezoid-weighted sum on a non-uniform one -- see `_period_dft`
        _dft = lambda B: self._period_dft(pss, B)
        if getattr(model, 'flicker', None) is None:
            groups = [lambda p: (self._sqrt_harmonics_of(model(wband(p)), _dft)
                                 if model is not None else
                                 self._cy_sqrt_harmonics(pss, wband(p)))]
        else:
            Pw = _dft(model.white)
            for l in ls:
                for lp in ls:
                    total += complex(rows[l] @ Pw[(lp - l) % Nn] @ np.conj(rows[lp]))
            groups = []
            for _key, Bc, EF in model.flicker:
                ef = self._uniform_exponent(Bc, EF)
                if ef is not None:
                    ## the element's SIGNED amplitudes where it states them
                    ## (any factor with `W W^dagger = B` serves the pair sum;
                    ## only this one knows the sign), else the PSD's root
                    _W = _signed.get(_key)
                    groups.append(lambda p, SB=(_dft(_W) if _W is not None else
                                                self._sqrt_harmonics_of(Bc, _dft)), ef=ef:
                                  (model.w1 / wband(p)) ** (0.5 * ef) * SB)
                else:
                    groups.append(lambda p, Bc=Bc, EF=EF: self._sqrt_harmonics_of(
                        Bc * (model.w1 / wband(p)) ** EF, _dft))
            for key in model.perband:
                groups.append(lambda p, key=key: self._sqrt_harmonics_of(
                    self._element_cy_samples(pss, wband(p))[key], _dft))
        for sqrt_at in groups:
            cache = {}
            ## every band the sum reaches, stacked once: BB[pi, k] =
            ## B_k^{(p)} with pi = p - pmin.  The (l, l') pair sum is then two
            ## fancy indexings and one einsum instead of N small products in
            ## Python.
            def _B(p, cache=cache, sqrt_at=sqrt_at):
                key = round(abs(f - p * f0) / f0, 12)
                if key not in cache:
                    cache[key] = sqrt_at(p)
                return cache[key]
            BB = np.asarray([_B(p) for p in range(pmin, pmax + 1)], dtype=complex)
            total = self._band_resolved_pairs(rows, ls, ks, Nn, pmin, BB, total)
        return float(np.real(total))

    @staticmethod
    def _band_resolved_pairs(rows, ls, ks, Nn, pmin, BB, total):
        for l in ls:
            for lp in ls:
                ## (B B^H)_j = sum_k B_k B_{k-j}^H: the partner index is
                ## k + l - l', NOT k + l' - l -- the mirror is invisible to
                ## the constant-modulation reduction (only k = 0 there).
                ## ⚠ NO CIRCULAR WRAP HERE: a partner beyond N/2 would be
                ## paired with the wrong BAND (each band carries its own
                ## weight), harmless in the white P-form and wrong here.
                ## History: `doc/shooting_history.md`,
                ## `PAC._band_resolved_pairs`.
                kp = ks + l - lp
                ok = np.abs(kp) <= Nn // 2
                kk, kk2 = ks[ok], kp[ok]
                pi = (l + kk) - pmin
                X = BB[pi, kk % Nn]
                Y = BB[pi, kk2 % Nn]
                Q = np.einsum('kij,klj->il', X, Y.conj())
                total += complex(rows[l] @ Q @ np.conj(rows[lp]))
        return float(np.real(total))

    def _cy_cycle_averaged(self, pss, w):
        """`CY` time-averaged over the orbit — Hull & Meyer's construction.

        ⚠⚠ VALID FOR GENTLE MODULATION ONLY, AND IT FAILS AS A FACTOR, NOT A
        PERCENTAGE: on a series switch + shunt capacitor it reads 16x a
        reference simulator at `goff/gon` = 1e-6 (1.000 unmodulated).  The
        averaged source injects `4kT <g>` for the WHOLE period, including
        the hold phase, where Hull & Meyer's own condition fails by six
        orders.  So this route is for a mixer's `gm` or a bias-dependent
        shot noise -- not for a switch.  A switch's noise is reachable
        exactly: `covariance` and `oscillator_covariance` evaluate `CY` at
        every step and need no averaging.

        This is what `_cy_reduced` refuses, done instead of refused, after
        Hull & Meyer (1993): *"cyclostationary noise sources, such as shot
        noise, may be modeled as MODULATED STATIONARY NOISE SOURCES.  The
        impulse response that is calculated INCLUDES THE EFFECT OF THIS
        MODULATION"*, with the stationary source at the cycle-averaged
        value.  The modulation is carried by the RESPONSE, `H_l`, rather
        than by the sources: ONE source per device.

        ⚠ Its condition: *"NONE OF THE LARGE-SIGNAL STATE VARIABLES MAY
        CHANGE SIGNIFICANTLY OVER THE DECAY TIME OF THE IMPULSE
        RESPONSE."*  A ringing (high-Q) impulse response breaks it, so it
        degrades as `lambda_2 -> 1`; `info` reports `|lambda_2|` so the
        caller can see which regime they are in.

        ⚠ SAMPLED ON THE ORBIT, NOT AT THE OPERATING POINT: `CY` is
        evaluated at every stored state and averaged with the step weights
        -- the same quadrature `diffusion_constant` uses, so the two remain
        comparable.

        History: `doc/shooting_history.md`, `PAC._cy_cycle_averaged`.
        """
        irn = pss.irefnode
        fp = pss.factored_period()
        tms = np.asarray(fp.times, dtype=float)
        T = float(fp.T)
        xs = np.asarray(pss.waveform[1], dtype=float)
        m = pss.cir.n - 1
        nsamp = min(len(tms) - 1, xs.shape[1])
        hs = self._period_weights(tms, nsamp, T, pss)
        acc = None
        for k in range(nsamp):
            xr = xs[:m, k]
            xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            cyk = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cyk,) = remove_row_col((cyk,), irn, pss.toolkit)
            cyk = np.asarray(cyk, dtype=complex) * hs[k]
            acc = cyk if acc is None else acc + cyk
        if acc is None:
            raise NotImplementedError(
                'PAC: the orbit has no stored samples to average CY over.')
        return acc / float(hs[:nsamp].sum())

    def _cy_reduced(self, pss, w):
        """`CY` with the reference node removed, refusing a moving one.

        ⚠ THE CHECK IS THE POINT.  A bias-dependent `CY` makes the sources
        CYCLOSTATIONARY, and then sidebands stop adding in power -- the
        stationary sum this class computes would be the wrong model, not
        merely an inaccurate one, and nothing downstream would say so.
        Sampled at three states on the converged orbit rather than argued
        from the element types, because a compact model's `CY` reads `x`
        and the discrete library's does not.

        ⚠⚠ IN EFFECT THIS REFUSES MOS pnoise: no physically correct MOS
        noise model has a state-independent `CY` (thermal `4kT gamma g_d0`
        with `g_d0` bias-dependent, flicker `I_D^AF`, gate shot `2qI_G`,
        trap rates reading the terminal voltages -- Mahmutoglu & Demir
        2015).  `PspMosLongChannel` is noiseless only by default (`fnt =
        0`).  The routes past the refusal are `cyclostationary=True` (exact
        for thermal and shot noise; flicker is the `|m|` fold, see
        `pnoise`) and `modulated=True` (Hull & Meyer's cycle average, see
        `_cy_cycle_averaged`).

        ⚠ This same check keeps the Ito/Stratonovich choice out of reach
        (`CY = GG^T`, so a state-dependent `CY` is a state-dependent `G`).
        Past it the two interpretations diverge, and Demir's escape -- "the
        noise signals are small compared with the deterministic signals" --
        may NOT carry for trap noise (a two-state Markov chain, not a small
        perturbation).  The tell would be a discrepancy in a MEAN but not in
        a variance.

        History: `doc/shooting_history.md`, `PAC._cy_reduced`.
        """
        irn = pss.irefnode
        fp = pss.factored_period()
        ## ⚠ THREE STATES ON THE ORBIT: the stored state half a period in,
        ## `x_last` and `x_prev`.  Never the zero vector -- it is on the
        ## orbit only by accident, and a switch model reading `goff` at
        ## v(ck) = 0 would refuse a linear time-invariant RC held by a DC
        ## clock as cyclostationary.
        _W = np.asarray(pss.waveform[1], dtype=float)
        _mid = np.delete(_W[:, _W.shape[1] // 2], irn, axis=0)
        states = [np.asarray(_mid, dtype=float).ravel()[:pss.cir.n - 1],
                  np.asarray(fp.x_last, dtype=float).ravel(),
                  np.asarray(fp.x_prev, dtype=float).ravel()[:pss.cir.n - 1]]
        mats = []
        for xr in states:
            xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
            cy = np.asarray(pss.cir.CY(xf, w), dtype=complex)
            (cy,) = remove_row_col((cy,), irn, pss.toolkit)
            mats.append(np.asarray(cy, dtype=complex))
        scale = max(float(np.max(np.abs(mats[0]))), 1e-300)
        for other in mats[1:]:
            drift = float(np.max(np.abs(other - mats[0]))) / scale
            if drift > 1e-9:
                raise NotImplementedError(
                    'PAC: this circuit has a BIAS-DEPENDENT CY (varies by '
                    '%.3g over the orbit), so its noise sources are '
                    'cyclostationary. The sidebands are then correlated '
                    'through the window Fourier coefficients and no longer '
                    'add in power, so the stationary sum here would be the '
                    'wrong model rather than an imprecise one. '
                    '⚠ THIS IS THE NORMAL CASE FOR ANY COMPACT MOS MODEL, '
                    'not an exotic one: thermal channel noise is '
                    '4kT.gamma.gd0 with gd0 bias-dependent, flicker goes '
                    'as I_D^AF, gate shot noise as 2qI_G, and trap capture '
                    'and emission rates read the terminal voltages. So '
                    'this is in effect a refusal of MOS pnoise, and the '
                    'route out is the CYCLOSTATIONARY construction rather '
                    'than a different device model -- BUILT: pass '
                    'cyclostationary=True (the PSD\'s own harmonics, exact; '
                    'modulated=True is the cycle average, which drops the '
                    'sideband correlation and read 0.53 of the truth on a '
                    'driven multiplier). Hull & Meyer (1993) '
                    'make it affordable -- ONE stationary source per '
                    'device at the CYCLE-AVERAGED current, with the '
                    'modulation carried by the impulse response, valid '
                    'while no large-signal state variable changes much '
                    'over the impulse response decay time. ⚠ THAT ROUTE '
                    'IS BUILT: pass modulated=True to use it. It is a '
                    'MODEL CHOICE with the validity condition above, not '
                    'a tolerance relaxation, which is why it is opt-in '
                    'and why this refusal is the default.'
                    % drift)
        return mats[0]

    ## How close to a harmonic of `f0` counts as "on" it, as a fraction of
    ## `f0` (`HARMONIC_GUARD`, below).  The deflated solve's conditioning is
    ## FLAT down to 1e-9 of `f0`, so the guard excludes only what has no
    ## finite answer: at an EXACT harmonic `1/(1 - alpha)` is a division by
    ## zero and the response is genuinely unbounded.
    ## History: `doc/shooting_history.md`, `PAC.HARMONIC_GUARD`.
    def _cy_at(self, pss, w, xr):
        """`CY` at ONE reduced state `xr`, reference row/column removed, with
        NO cyclostationarity check -- for the routes that evaluate the
        source at every step and therefore model a modulated source
        exactly (`covariance`, `oscillator_covariance`).  The stationary
        sum in `pnoise` cannot, which is why `_cy_reduced` refuses there.
        """
        irn = pss.irefnode
        xr = np.asarray(xr, dtype=float).ravel()[:pss.cir.n - 1]
        xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
        cy = np.asarray(pss.cir.CY(xf, w), dtype=complex)
        (cy,) = remove_row_col((cy,), irn, pss.toolkit)
        return np.asarray(cy, dtype=complex)

    def _lyap_cy(self, pss, w, xr):
        """`CY` as the Lyapunov pieces read it: `_cy_at`, except while a
        COLOURED covariance runs (`_white_cy` set by `covariance` /
        `event_jitter`), when it is the WHITE part `A(x)` of the component
        model alone -- the coloured part is added in the frequency domain
        (`_coloured_covariance`), and `CY` at `w0` would count it again."""
        white = getattr(self, '_white_cy', None)
        if white is None:
            return self._cy_at(pss, w, xr)
        return white(xr)

    HARMONIC_GUARD = 1e-12

    def _check_circuit(self, pss):
        """The operating point must belong to THIS circuit, not a similar one.

        ⚠ A DRIVEN OSCILLATOR MAKES THIS A CORRECTNESS TRAP RATHER THAN A
        TYPO GUARD.  The natural way to model one is to solve the PSS of
        the bare oscillator and then treat the injection as a perturbation
        -- and it is wrong, because the injection DEVICE is present even
        when its SIGNAL is zero.  Buonomo & Lo Schiavo: "in absence of the
        injection signal, the injection circuit affects the basic LC
        oscillator by CHANGING THE NONLINEARITY OF THE FEEDBACK LOOP ...
        [it] can affect the start-up condition of the basic differential LC
        oscillator OR ITS OSCILLATION AMPLITUDE, or both."

        So the free-running orbit of the circuit-with-the-device is not the
        orbit of the circuit-without-it, and every Floquet quantity built
        on the wrong one inherits the error -- monodromy, PPV, phase noise.
        The analysis would converge and report a plausible number.

        The reference-node check below catches a mismatched `refnode`; it
        cannot catch this, because two circuits differing by one device
        have the same reference node and often the same node count.
        """
        if self.cir is not pss.cir:
            raise ValueError(
                'PAC: this analysis was built on a different circuit object '
                'than the PSS it was handed. If that is deliberate -- e.g. '
                'solving the PSS of a bare oscillator and perturbing a '
                'version with an injection device added -- it is a '
                'CORRECTNESS error, not a bookkeeping one: the injection '
                'device changes the free-running orbit even with its signal '
                'at zero, so the base solution is the wrong one to '
                'linearise about. Solve the PSS on the SAME circuit.')

    def _check_harmonic(self, pss, freq, what):
        """Refuse an autonomous small-signal solve sitting on a harmonic.

        ⚠ `I - exp(-j w T) M` IS SINGULAR AT EVERY HARMONIC OF `f0`, NOT
        JUST AT DC, AND ONLY FOR AN OSCILLATOR.  At `w = k w0` the factor
        `exp(-j w T)` is 1 and the operator is `I - M`, which an autonomous
        circuit's unit multiplier makes singular; `sigma_min` falls LINEARLY
        with the distance to the nearest harmonic.  A driven circuit has no
        unit multiplier and no singularity, harmonics included.

        ⚠ AND IT IS PHYSICS, NOT CONDITIONING.  A perturbation at a
        harmonic is a perturbation along the orbit, and an oscillator's
        response to that is unbounded phase drift -- there is no bounded
        periodic answer to return.  So this refuses rather than tightening
        a tolerance, and says which quantity to ask for instead.

        History: `doc/shooting_history.md`, `PAC._check_harmonic`.
        """
        if not getattr(pss, 'autonomous', False):
            return
        f0 = 1.0 / float(pss.factored_period().T)
        r = abs(float(freq)) / f0
        d = abs(r - round(r))
        if d <= self.HARMONIC_GUARD:
            raise ValueError(
                'PAC: %s at %.6g Hz is on harmonic %d of this OSCILLATOR\'s '
                'own frequency (%.6g Hz), where I - exp(-j w T) M is '
                'singular -- the unit Floquet multiplier makes it exactly '
                'I - M there. That is physics, not conditioning: a '
                'perturbation along the orbit produces unbounded phase '
                'drift, so there is no bounded periodic response to '
                'return. Ask off-harmonic, or ask for the phase quantity '
                'instead (PSS.ppv()).' % (what, float(freq), round(r), f0))

    @staticmethod
    def _gmres_checked(A, b, rtol, what):
        """GMRES (`_arnoldi_gmres`), judged by its RESIDUAL.

        These operators are `2m x 2m` and often tiny, so the Krylov space is
        exhausted in a handful of steps; the next vector is then numerically
        zero, which is a LUCKY breakdown -- the solution is exact -- and a
        status flag (scipy's `info = 4`) reports it as a failure.  So the
        residual decides, and it is also the tolerance decision.  A genuine
        failure still fails, with the residual quoted, because the real
        cause near a harmonic is that the operator is nearly singular there
        and no tolerance will fix it.

        History: `doc/shooting_history.md`, `PAC._gmres_checked`.
        """
        n = b.shape[0]
        x, relres, _H, _k = _arnoldi_gmres(
            A.matvec, b, rtol=rtol, maxiter=min(n, 200))
        info = 0
        r = float(np.linalg.norm(b - A.matvec(x)))
        scale = max(float(np.linalg.norm(b)), 1e-300)
        if r / scale > max(1e3 * rtol, 1e-8):
            raise RuntimeError(
                'PAC: %s did not converge (info=%r, relative residual '
                '%.3e). Near a harmonic of the oscillator this operator is '
                'genuinely near-singular and a smaller tolerance will not '
                'help -- move the offset, or ask for the phase quantity.'
                % (what, info, r / scale))
        return x

    def _refuse_coloured(self, pss, what):
        """Refuse a coloured source where the machinery assumes WHITE.
        (`covariance` and `event_jitter` take a band instead and do not
        come here with one -- `_coloured_prepare`.)

        ⚠ THE TRAP IS THAT NOTHING ELSE WOULD OBJECT. The Lyapunov
        recursion, `diffusion_constant` and eq (22)'s collapse all read
        `CY` at ONE frequency and treat it as the noise intensity at every
        frequency; a coloured source folded that way returns a plausible
        number, not an error (A4d names exactly this shape). Detected by
        evaluating the reduced `CY` at two frequencies -- colour is
        frequency dependence, bias dependence is what `_cy_reduced`
        refuses separately.
        """
        if self._coloured_present(pss):
            raise NotImplementedError(
                'PAC.%s: a noise source in this circuit is COLOURED (its CY '
                'differs between w0 and 10 w0), and this routine assumes '
                'white sources -- it would fold CY at one frequency as if '
                'it held at every frequency and return a plausible wrong '
                'number. Use the frequency-resolved surfaces (pnoise, '
                'sampled_variance, phase_psd/coloured_diffusion on an '
                'oscillator, where 1/f noise makes the phase growth '
                'non-diffusive), covariance/event_jitter with a band '
                '(fmin, fmax) on a driven circuit, or the white-through-filter '
                'form of the source.' % what)

    def _coloured_present(self, pss):
        """Whether any noise source of the circuit is COLOURED -- see
        `_refuse_coloured`, which asks this."""
        w1 = 2.0 * np.pi / float(pss.period)
        ## ⚠ Colour is asked at fixed state, two frequencies: it is separable
        ## from the bias question, and asking it through `_cy_reduced` would
        ## refuse every MODULATED source the covariance routes (which
        ## evaluate `CY` per step) handle exactly.
        ## ⚠⚠ PER ENTRY, AND AT MORE THAN ONE STATE: each entry is judged
        ## against ITS OWN magnitude (a drain's flicker colour sits orders
        ## below a gate resistor's white 4kT/rg in the same matrix), at
        ## `x_last` and at states spread over the orbit (a switch OFF at
        ## t = 0 is white there and coloured elsewhere).  Exact zeros and
        ## white entries give identical values at both frequencies, so
        ## neither can fire.
        ## `_cy_at` takes a REDUCED state: the reference row is removed here.
        ## History: `doc/shooting_history.md`, `PAC._refuse_coloured`.
        _xl = np.asarray(pss.factored_period().x_last, dtype=float).ravel()
        _states = [_xl]
        _wf = getattr(pss, 'waveform', None)
        if _wf is not None:
            _W = np.delete(np.asarray(_wf[1], dtype=float), pss.irefnode,
                           axis=0)
            for _k in sorted(set(np.linspace(0, _W.shape[1] - 1,
                                             8).astype(int))):
                _states.append(_W[:, _k])
        for _xr in _states:
            c1 = self._cy_at(pss, w1, _xr)
            c2 = self._cy_at(pss, 10.0 * w1, _xr)
            den = np.maximum(np.abs(c1), np.abs(c2))
            if np.any(np.abs(c1 - c2) > 1e-9 * den):
                return True
        return False

    def _lyapunov_pieces(self, pss, what):
        """The per-step maps, injections and one-period accumulation.

        Returns `(As, Qs, K1, M, m, n)`: the step maps `A_j`, the noise
        injections `Q_j`, the covariance `K1` reached after one period
        starting from zero, the monodromy `M`, and the two widths.

        ⚠ SHARED BY THE DRIVEN AND AUTONOMOUS ROUTES ON PURPOSE.  The two
        differ only in what they do with `I - M kron M`: `covariance`
        inverts it, `oscillator_covariance` borders it because it is
        singular there.  Everything upstream -- the `CY/2` convention, the
        `b = 0` restriction, the `C` ring the forward recursion sees -- is
        one implementation, so the pair cannot drift apart over exactly this
        factor of two.

        History: `doc/shooting_history.md`, `PAC._lyapunov_pieces`.
        """
        if getattr(self, '_white_cy', None) is None:
            ## (a coloured `covariance` / `event_jitter` has set the WHITE
            ## part for the pieces and adds the coloured one itself)
            self._refuse_coloured(pss, what)
        fp = pss.factored_period()
        if fp.is_glm:
            ## reached with monodromy='native' only (`_lyapunov_host`)
            return self._lyapunov_pieces_glm(pss, fp.state_map(), what)
        if fp.is_stage:
            ## the stage method's per-step map + its stage injection (or the
            ## SAME exact Van Loan integral) -- see `_lyapunov_pieces_stage`
            return self._lyapunov_pieces_stage(pss, fp, what)
        if not fp.is_pair:
            return self._lyapunov_pieces_plain(pss, fp, what)
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
                                     _W[:, min(k + 1, _W.shape[1] - 1)]))
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

    def _lyapunov_pieces_plain(self, pss, fp, what):
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
                                     _W[:, min(k + 1, _W.shape[1] - 1)]))
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

    def _lyapunov_pieces_glm(self, pss, sm, what):
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
                        pss, w0, np.delete(np.asarray(y, dtype=float), irn)),
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

    def _lyapunov_pieces_stage(self, pss, fp, what):
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
        if getattr(self, '_white_cy', None) is None:
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
            CYn = self._lyap_cy(pss, w0, xk)
            A_k = np.column_stack([
                np.asarray(pss._monodromy_matvec_stage([step], e), dtype=float)
                for e in np.eye(m)])
            As.append(A_k)
            Q_k = self._stage_injection(pss, fp, k, w0)
            Qs.append(Q_k if Q_k is not None
                      else self._vanloan_step_injection(Cn, Gn, CYn, hs[k]))
        K = np.zeros((n, n))
        for A_k, Q_k in zip(As, Qs):
            K = A_k @ K @ A_k.T + Q_k
        M = np.column_stack([np.asarray(fp.matvec(e), dtype=float)
                             for e in np.eye(n)])
        return As, Qs, K, M, m, n

    def _stage_injection(self, pss, fp, k, w):
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
                CYi = np.real(np.asarray(self._lyap_cy(pss, w, yi), dtype=complex))
                Q += wts[i] * self._vanloan_step_injection(
                    np.asarray(pss._C_at(yi), dtype=float),
                    np.asarray(pss._G_at(yi), dtype=float), CYi, h)
            return 0.5 * (Q + Q.T)
        Q = np.zeros((m, m))
        for i in range(s):
            yi = np.delete(np.asarray(states[k * s + i], dtype=float), irn)
            CYi = np.real(np.asarray(self._lyap_cy(pss, w, yi), dtype=complex))
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

    def _orbit_rate(self, pss, event_nodes):
        """`xdot` at every node of the solved orbit (reduced width) -- THE
        DAE'S OWN DERIVATIVE at the node's state: on the
        differential rows ``C(x) xdot = -(i(x) + u(t))``, on the algebraic
        rows (a zero row of `C`) the differentiated constraint ``G(x) xdot
        = -du/dt``; one small solve per node, exact for the discrete state
        and independent of the step.  The rate converts a node's motion in
        time into a state change: see `_fixed_time_event_columns`.

        The three-node stencil (`_orbit_rate_stencil`: second order, and
        one-sided at a landed event, where a step cannot fit a parabola to
        a fast exponential) is only the fallback where the assembled matrix
        is singular (an index above one), and says so.

        History: `doc/shooting_history.md`, `PAC._orbit_rate`."""
        ts = np.asarray(pss.waveform[0], dtype=float)
        X = np.delete(np.asarray(pss.waveform[1], dtype=float),
                      pss.irefnode, axis=0)
        N = len(ts) - 1
        m = X.shape[0]
        out = np.zeros((N + 1, m))
        analysis = getattr(pss.par, 'analysis', None)
        ok = True
        for j in range(N + 1):
            ## ⚠ NODE 0 IS EVALUATED AS NODE N.  A source's derivative at
            ## exactly its start (`VSin` clamps `t - td` at 0: SPICE's rule
            ## for a transient) is the LEFT one -- 0 -- where the periodic
            ## steady state, t = 0 == t = T, has the right one.  Node N is
            ## the same point.
            jj = N if j == 0 else j
            x = X[:, jj]
            t = float(ts[jj])
            try:
                C = np.asarray(pss._C_at(x), dtype=float)
                k = np.asarray(pss._k_at(x, t), dtype=float).ravel()
                alg = [i for i in range(m) if not np.any(C[i, :])]
                A = C.copy()
                b = k.copy()
                if alg:
                    G = np.asarray(pss._G_at(x), dtype=float)
                    ud = np.delete(np.asarray(pss.cir.dudt(t, analysis=analysis),
                                              dtype=float).ravel(), pss.irefnode)
                    A[alg, :] = G[alg, :]
                    b[alg] = -ud[alg]
                out[j] = np.linalg.solve(A, b)
            except (np.linalg.LinAlgError, ValueError):
                ok = False
                break
        if ok:
            return out
        warnings.warn(
            'PAC._orbit_rate: the DAE derivative could not be assembled at a '
            'node (a singular differential/algebraic split -- an index above '
            'one?); falling back to the three-node stencil, which is second '
            'order in the step and one-sided at a landed event.',
            RuntimeWarning, stacklevel=3)
        return self._orbit_rate_stencil(pss, event_nodes)

    def _orbit_rate_stencil(self, pss, event_nodes):
        """The three-node parabola: one-sided AT a landed event and at the
        node after one, central elsewhere, periodic at the ends.  Kept as
        `_orbit_rate`'s fallback.

        History: `doc/shooting_history.md`, `PAC._orbit_rate_stencil`."""
        ts = np.asarray(pss.waveform[0], dtype=float)
        X = np.delete(np.asarray(pss.waveform[1], dtype=float),
                      pss.irefnode, axis=0)
        N = len(ts) - 1
        T = float(ts[-1] - ts[0])
        ev = set(int(j) for j in event_nodes)
        out = np.zeros((N + 1, X.shape[0]))

        def _t(i):
            k = i % N
            return float(ts[k]) + T * ((i - k) // N)

        for j in range(N + 1):
            if j in ev or (j % N) in ev:
                a, b, c = j - 2, j - 1, j
            elif (j - 1) in ev or ((j - 1) % N) in ev:
                a, b, c = j, j + 1, j + 2
            else:
                a, b, c = j - 1, j, j + 1
            ta, tb, tc = _t(a), _t(b), _t(c)
            xa, xb, xc = X[:, a % N], X[:, b % N], X[:, c % N]
            tj = _t(j)
            out[j] = (xa * ((tj - tb) + (tj - tc)) / ((ta - tb) * (ta - tc))
                      + xb * ((tj - ta) + (tj - tc)) / ((tb - ta) * (tb - tc))
                      + xc * ((tj - ta) + (tj - tb)) / ((tc - ta) * (tc - tb)))
        return out

    def _fixed_time_event_columns(self, pss):
        """The event columns at every node AT FIXED TIME: ``Pk_j - xdot_j
        tau_j^T`` with `tau_j = sum_{i<j} dh_i/dtheta` (seconds) the node's
        own motion when the crossings move (`_event_remap` scales the steps
        between two crossings together) and `xdot_j` from `_orbit_rate`.
        ⚠ The stored `Pk_nodes` are the derivatives of "the state at node
        j", a point whose TIME moves with theta; a consumer that reports
        a response or a covariance at the grid's times needs this form --
        without it a noiseless ramp source reads 0.225 kT/C of a
        threshold's noise (its node's share of the moving segment) and the
        sideband response along a staged oscillator's orbit is O(1) off
        the exact one while its period node is exact.  Returns
        ``(Pk_fixed, tau, xdot)``."""
        ev = pss._event_columns
        th = np.asarray(pss._state_event_fracs, dtype=float)
        _fr, hsens, _nd = pss._event_remap(
            np.asarray(pss._grid_fracs, dtype=float), th, th, float(pss.period))
        tau = np.vstack((np.zeros((1, len(th))),
                         np.cumsum(np.asarray(hsens, dtype=float), axis=0)))
        xdot = self._orbit_rate(pss, ev['nodes'])
        Pk = np.asarray(ev['Pk_nodes'], dtype=float)
        return Pk - xdot[:, :, None] * tau[:, None, :], tau, xdot

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

        Returns ``{'sigma': (K,) s, 'cov_fraction': (K, K), 'fractions':
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
        try:
            As, Qs, K1, M, m, n = self._lyapunov_pieces(pss, 'covariance')
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
        finally:
            self._white_cy = None
        if col is not None:
            ## the crossings' coloured motion, from the same bordered forced
            ## responses (`_forced_responses`' shifts)
            _Kc, Dc = self._coloured_covariance(pss, col, m, n)
            cov = cov + Dc
        cov = 0.5 * (cov + cov.T)
        T = float(pss.period)
        return {'sigma': np.sqrt(np.clip(np.diag(cov), 0.0, None)) * T,
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
        coloured components and the band, and the Lyapunov pieces are set to
        read the WHITE part of each source (`_white_cy`, `_lyap_cy`) -- the
        caller clears it.  Refuses what cannot be integrated: no `fmin` (a
        1/f variance grows as ``ln(fmax/fmin)`` without limit), a component
        that is not a power law, a circuit whose `CY` is not the sum of its
        elements'."""
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
        model = self._cy_components_model(pss, fmin, f0, states)
        if model is None:
            raise NotImplementedError(
                'PAC.%s: this circuit\'s CY is not the sum of its elements\', '
                'so its coloured part cannot be separated from the white one '
                '(see the warning above).' % what)
        if model.perband:
            raise NotImplementedError(
                'PAC.%s: the noise of %s is coloured but not a power law '
                '(thermal plus 1/f^EF), so there is no density to integrate '
                'over the band. Use sampled_noise / pnoise, which evaluate '
                'it per band.' % (what, ', '.join('.'.join(k) for k in model.perband)))
        self._warn_signed_unused(model, 'PAC.%s' % what)
        amp = getattr(model, 'amplitude', None) or {}
        comps = []
        for key, B, EF in model.flicker:
            ef = self._uniform_exponent(B, EF)
            if ef is None:
                raise NotImplementedError(
                    'PAC.%s: the coloured noise of %s carries different '
                    'power-law exponents in different entries, so it has no '
                    'one amplitude to replay. Use sampled_noise / pnoise.'
                    % (what, '.'.join(key)))
            ## ⚠ THE SIGN: the element's stated amplitudes where it has them
            ## (`W W^H = B` with the sign of the modulation); `sqrt(B)` is
            ## the sign-blind |m| process (`_warn_signed_unused` said so)
            W = amp.get(key)
            W = np.asarray(W if W is not None else self._psd_sqrt(B), dtype=complex)
            comps.append((key, W, float(ef)))
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
                    mdl = self._cy_components_model(pss, fmin, f0, states=[xr])
                cache[key] = np.asarray(mdl.white[0], dtype=complex)
            return cache[key]
        self._white_cy = white
        return {'fp': fp, 'counts': counts, 'comps': comps, 'w1': model.w1,
                'fmin': fmin, 'fmax': fmax, 'ppd': int(points_per_decade)}

    def _coloured_covariance(self, pss, col, m, n):
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
        log grid, `points_per_decade` as `sampled_variance`'s, and the
        trapezoid in `ln(nu)`.  Gear's PAIR covariance carries `(x_j,
        x_{j-1})`: the previous node's response, node -1 being node N - 1
        a period back (``e^{-j 2 pi nu T}``).

        ⚠ A SLOPE, NOT A STATE: no shaping filter, no fitted Lorentzian
        ladder -- the exact power law over a hard band, as `sampled_variance`
        and `coloured_diffusion` take it.

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
        ## the trapezoid in ln(nu), not in nu: `int F dnu = int nu F dln(nu)`
        ## is EXACT for a pure 1/f density and spectrally accurate for one
        ## whose response flattens at both band ends (the linear trapezoid
        ## on the same log grid carries (r - 1)^3 / 6 per point, r the
        ## grid ratio: 6e-4 of a 1/f band at 40 per decade)
        nn = max(2, int(np.ceil(col['ppd'] * np.log10(fmax / fmin)))) + 1
        nus = np.geomspace(fmin, fmax, nn)
        dl = np.diff(np.log(nus))
        wq = np.zeros(nn)
        wq[:-1] += 0.5 * dl
        wq[1:] += 0.5 * dl
        wq = wq * nus
        offs = np.concatenate(([0], np.cumsum(col['counts'])))
        K = np.zeros((N + 1, n, n))
        D = None
        zero = np.zeros(m, dtype=complex)
        for _key, W, ef in col['comps']:
            for s_ in range(W.shape[2]):
                Wc = W[:, :, s_]
                if not np.any(Wc):
                    continue
                u_points = [Wc[offs[j]:offs[j + 1]] for j in range(N)]
                ys, dths = self._forced_responses(pss, fp, nus, zero,
                                                  u_points=u_points)
                for i, nu in enumerate(nus):
                    dens = wq[i] * (col['w1'] / (2.0 * np.pi * nu)) ** ef
                    y = np.asarray(ys[i], dtype=complex)[:N + 1]
                    if n != m:
                        prev = np.vstack((y[N - 1:N] * np.exp(-2j * np.pi * nu * T),
                                          y[:N]))
                        y = np.hstack((y, prev))
                    K += dens * np.real(np.einsum('ji,jk->jik', y, y.conj()))
                    if dths[i] is not None:
                        dd = np.asarray(dths[i], dtype=complex)
                        Di = dens * np.real(np.outer(dd, dd.conj()))
                        D = Di if D is None else D + Di
        K = 0.5 * (K + np.swapaxes(K, 1, 2))
        return K, D

    def covariance(self, pss, samples=False, fmin=None, fmax=None,
                   points_per_decade=40):
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

        Returns `K0`, the covariance at `t = 0`; with `samples=True`,
        `(K0, [K_j])`, the covariance at every step, which is the
        time-varying statistic this exists to produce.

        ⚠ A COLOURED SOURCE NEEDS A BAND (2026-09-25).  With a 1/f source
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
        Refused: a colour that is not a power law, a component whose
        exponent differs between its entries, and trap's plain map (its
        covariance is on the (x, iq) pair).

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
        try:
            As, Qs, K1, M, m, n = self._lyapunov_pieces(pss, 'covariance')
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
        finally:
            self._white_cy = None
        if col is not None:
            Kc, _dth = self._coloured_covariance(pss, col, m, n)
            K0 = K0 + Kc[0]
            if seq is not None:
                seq = [a + b for a, b in zip(seq, Kc)]
        if not samples:
            return K0
        return K0, seq

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
        white part and power-law part separately; `_cy_components_model`).
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

        ⚠ POWER LAW BETWEEN THE POINTS (`_loglog_integral`, 2026-09-25): the
        density is interpolated linearly in log-log and each interval
        integrated exactly, so a white band and any pure power law are
        exact and the error is second order where the density BENDS.  The
        linear trapezoid it replaced overestimates a 1/f band by
        `(r - 1)^3 / 6` per point (r the grid ratio); on the switched
        sampler at 40 per decade, against a 640-per-decade reference:
        1/f +3.7e-4 / +5.1e-4 before, +6.6e-5 / +6.1e-5 now; white +
        1/f +7.7e-5 / +5.1e-4 before, +2.6e-5 / +6.1e-5 now; white 1e-11
        both.  (The trapezoid in ln f, exact for 1/f, is +2.8e-4 on a
        white band -- not used.)

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

    @staticmethod
    def _loglog_integral(S, fs):
        """``int S df`` over the grid `fs` (last axis of `S`), `S` a power
        law between neighbouring points: on `[f1, f2]`, ``S = S1
        (f/f1)^p`` through both ends, integrated exactly --
        ``S1 f1 ln(r) (e^z - 1)/z``, ``z = ln(S2 f2 / (S1 f1))``, stable at
        ``z -> 0`` (1/f).  An interval with an end that is not positive
        (no power law passes through it) takes the linear trapezoid.

        History: `doc/shooting_history.md`, `PAC.sampled_variance`."""
        S = np.asarray(S, dtype=float)
        f1, f2 = fs[:-1], fs[1:]
        S1, S2 = S[..., :-1], S[..., 1:]
        L = np.log(f2 / f1)
        ok = (S1 > 0.0) & (S2 > 0.0)
        with np.errstate(divide='ignore', invalid='ignore'):
            z = np.log(np.where(ok, S2, 1.0) / np.where(ok, S1, 1.0)) + L
            phi = np.where(np.abs(z) < 1e-8, 1.0 + 0.5 * z,
                           np.expm1(z) / np.where(z == 0.0, 1.0, z))
        pw = S1 * f1 * L * phi
        lin = 0.5 * (S1 + S2) * (f2 - f1)
        return np.sum(np.where(ok, pw, lin), axis=-1)

    def jitter_metrics(self, pss, output, time, fmin, fmax, kmax=4,
                       maxsidebands=None, nfreq=601, dc_rectangle=False):
        """Edge jitter at ONE instant: `sigma_t`, the across-period
        correlation `rho_k`, and the three metrics that are functions of it.
        DRIVEN circuits (an oscillator is refused by `_sampled_series`; its
        accumulating share is `diffusion_constant`'s `c`).

        `sampled_noise` returns the one-sided PSD of the SAMPLE SERIES
        `y(t0 + kT)`, so that series' own autocovariance is its cosine
        transform -- no new machinery, and none of Demir 1996:

            R_k   = int_fmin^fmax S(f; t0) cos(2 pi f k T) df
            rho_k = R_k / R_0

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
        rect = float(S[0]) * fmin if dc_rectangle else 0.0
        R = np.array([float(trapezoid(S * np.cos(2.0 * np.pi * fs * k * T), fs))
                      + rect for k in range(kmax + 1)])
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
        model = self._cy_components_model(pss, float(np.min(fr)), f0, states)
        white, scaled, perband = [], [], []
        if model is None:
            ## the elements do not sum to the circuit's CY (warned): one
            ## joint component over the whole circuit
            perband.append(lambda w: self._cy_at_states(pss, w, states))
        else:
            white = [self._psd_sqrt(A) for _key, A in model.white_parts]
            self._warn_signed_unused(model, 'PAC.sampled_noise')
            for _key, Bc, EF in model.flicker:
                ef = self._uniform_exponent(Bc, EF)
                if ef is not None:
                    _W = (getattr(model, 'amplitude', None) or {}).get(_key)
                    scaled.append((_W if _W is not None else self._psd_sqrt(Bc), ef))
                else:
                    perband.append(lambda w, Bc=Bc, EF=EF: Bc * (model.w1 / w) ** EF)
            for key in model.perband:
                perband.append(lambda w, key=key: self._element_cy_samples(
                    pss, w, states)[key])

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
                for SB, ef in scaled:
                    R = E @ np.einsum('ji,jik->jk', Sv, SB)
                    c = (model.w1 / (2.0 * np.pi * np.abs(nu))) ** ef
                    _pb += c * np.sum(np.abs(R) ** 2, axis=1)
                for comp in perband:
                    for bi, nb in enumerate(nu):
                        R = E[bi] @ np.einsum(
                            'ji,jik->jk', Sv,
                            self._psd_sqrt(comp(2.0 * np.pi * abs(nb))))
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

    def _stage_times(self, pss, fp):
        """The stage abscissae `t_j + c_k h_j` of one period, step by step,
        `(N s,)` -- the injection times of a stage method's source.  A GLM's
        map on the state lists its steps' own (`_GLMStateStep`: the stages,
        then the substages of the startup that opens a step)."""
        tms = np.asarray(fp.times, dtype=float)
        if fp.is_glm:
            return np.asarray([t for j, st in enumerate(fp.step_objects())
                               for t in st.injection_times(tms[j])],
                              dtype=float)
        out = []
        for j, st in enumerate(fp.steps):
            h = tms[j + 1] - tms[j]
            out.extend(tms[j] + st.c * h)
        return np.asarray(out, dtype=float)

    def _stage_states(self, pss, fp):
        """The stage states of one period, `N s` full-width vectors in the
        order of `_stage_times` -- one re-traversal of the converged orbit
        on a FRESH inner transient (the run's own is restored), cached per
        factored period.  A GLM's steps carry theirs (the walk stored them)."""
        if fp.is_glm:
            return [y for st in fp.step_objects()
                    for y in st.injection_states()]
        cache = getattr(pss, '_sampled_stage_cache', None)
        if cache is not None and cache[0] is fp:
            return cache[1]
        tms = np.asarray(fp.times, dtype=float)
        hs = np.diff(tms)
        m = pss.cir.n - 1
        x = np.asarray(pss._period_state[1], dtype=float)[:m]
        saved = getattr(pss, '_tran', None)
        pss._tran = pss._new_transient(pss._integrator_for(pss.par.method))
        try:
            pss._begin_period(x)
            tr = pss._transient()
            Ys = []
            for j in range(len(fp.steps)):
                x = np.asarray(pss.solve_timestep(x, tms[j + 1], hs[j]),
                               dtype=float)
                Ys.extend(np.asarray(y, dtype=float) for y in tr._rk_Y)
        finally:
            pss._tran = saved
        pss._sampled_stage_cache = (fp, Ys)
        return Ys

    def _stage_pass(self, pss, fp, lam0, seed=None):
        """One reverse pass of a stage period map (`dirk` or `full`, or a
        GLM's map on the state): returns the final costate and the coupling
        vectors, one per injection point (`_stage_times`; `h sum_i A_ik p_i`
        per stage of a stage method) -- the sensitivity of the costate's
        functional to a unit source there is minus that (see
        `_sideband_forced`, whose loop this is).  `seed = (k0, v)` adds `v`
        to the costate on the state after step `k0`'s update: the output at
        `t_{k0}` couples to the sources of earlier steps only."""
        steps = fp.step_objects()
        N = len(steps)
        lam = fp.extract_T(np.asarray(lam0, dtype=complex).copy())
        per = [None] * N
        for j in range(N - 1, -1, -1):
            st = steps[j]
            lam, r = st.adjoint(lam)
            per[j] = st.couplings(r)
            if seed is not None:
                if isinstance(seed, dict):
                    if j in seed:
                        lam = fp.inject(lam, np.asarray(seed[j], dtype=complex))
                elif j == seed[0]:
                    lam = fp.inject(lam, seed[1])
        return (fp.seed_T(lam),
                np.asarray([v for row in per for v in row], dtype=complex))

    @staticmethod
    def _psd_sqrt(Cs):
        """Symmetric PSD square roots of a stack `(..., n, n)`."""
        Cs = np.asarray(Cs, dtype=complex)
        Cs = 0.5 * (Cs + np.conj(np.swapaxes(Cs, -1, -2)))
        lam, U = np.linalg.eigh(Cs)
        return np.einsum('...ik,...k,...jk->...ij', U,
                         np.sqrt(np.clip(np.real(lam), 0.0, None)), U.conj())

    ## ⚠ ON AN ASYMMETRIC ORBIT this route is validated by a direct SDE
    ## simulation of the variational system (the calibrated `Var(i) =
    ## CY/(2h)` injection, phase projected out every step), which shares no
    ## Lyapunov solve and no modal sum: 0.02 % at `a = 0.30` on van der Pol
    ## + `a u^2`.  ⚠ Asymmetry changes the MODE SHAPES, so the noise
    ## projected onto the orbital direction grows and the transverse
    ## variance RISES despite the faster relaxation (`|lam2|` falls) -- the
    ## physical argument that it should shrink is wrong.
    ## History: `doc/shooting_history.md`, `PAC.oscillator_covariance`.
    def oscillator_covariance(self, pss, samples=False):
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
        As, Qs, K1, M, m, n = self._lyapunov_pieces(
            pss, 'oscillator_covariance')
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
            ## ⚠ THE TRAPEZOIDAL PLAIN PAIR `(x, iq)` (2026-09-24; refused
            ## before): its map re-seeds `iq` at every period start, so its
            ## last `m` columns are zero and its null vectors follow from the
            ## state map's -- the left one `[v; 0]`, the right one the tangent
            ## with the `iq` block it carries, `M[:, :m] u`.  The result is as
            ## good as trap's own map: first order on a limit cycle (d -2.5 %
            ## at 200 points, -1.1 % at 800 on van der Pol, against a radau
            ## reference; the trbdf2 twin -4.8e-6).
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

        ⚠ `k_cycle` IS THE LARGE-`k` FORM, with `rho_k` taken to zero.  The
        orbital part's across-period correlation is not computed here, so at
        small `k` the true k-cycle jitter is LOWER than this returns (the
        `(1 - rho_k)` factor is below 1).  Treat small-`k` entries as an upper
        bound.  For a DRIVEN circuit use `jitter_metrics`, which computes
        `rho_k` properly from the sample series.

        ⚠ THE INSTANT IS THE CALLER'S.  This does not hunt for a crossing: a
        threshold taken from a simulated record can be biased by startup, and
        that moves the instant off the steepest point.  `instant` in the
        result is the grid point used.

        Returns a dict: `sigma_t`, `A` (= sigma_t^2), `c`, `slew`, `k_cycle`
        (k = 1..kmax), `instant`, `d`, `projection_share`.

        History: `doc/shooting_history.md`, `PAC.oscillator_edge_jitter`.
        """
        import warnings as _warnings
        self._check_circuit(pss)
        ## ONE ORBIT: the covariance's host (a GLM's or trap's twin) supplies
        ## the period, the factored period and the PPV as well -- read off
        ## the run itself they came from another discretisation (trap's own
        ## period against its twin's covariance, until 2026-09-24)
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
        j = int(np.argmin(np.abs(times - float(time))))

        ## ⚠ THE SLOPE IS TAKEN AT THE REQUESTED INSTANT, NOT AT THE SNAPPED
        ## GRID POINT, and the difference is the whole correction.  `time` is
        ## typically a threshold crossing, which sits at a different FRACTION
        ## of a step on every grid; differentiating at the nearest sample
        ## instead reproduces the straddling value, and every quantity here
        ## goes as 1/s^2.
        lo = max(min(j - 2, nt - 5), 0)
        tt = times[lo:lo + 5] - float(time)
        a2, b2, _c2 = np.polyfit(tt, row[lo:lo + 5], 2)
        slew = float(b2)
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
            'k_cycle': np.sqrt(c * ks * T + 2.0 * A),
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
        K_orb, _d, _info = self.oscillator_covariance(pss)
        K = np.asarray(K_orb, dtype=float)
        n = K.shape[0]
        modes = pss.floquet_modes(pss, nmodes=(n if nmodes is None
                                               else int(nmodes)))
        V = np.column_stack([m['v0'] for m in modes])
        cw = V.conj().T @ K @ V
        return cw, modes, K

    ORBITAL_HARMONICS = 32

    ## Half-wave asymmetry above which `orbital_correlation` warns of its O(h)
    ## residual and `orbital_spectrum` of its over-statement; 0.02 is a decade
    ## inside the smallest asymmetry at which the error was visible.
    ORBITAL_ASYMMETRY_LIMIT = 0.02

    def _orbit_asymmetry(self, pss):
        """Half-wave asymmetry of the orbit, in [0, ~1].

        `max|x(t) + x(t + T/2)| / max|x|` on the first state row -- zero for a
        half-wave symmetric orbit (van der Pol), growing as the orbit
        distorts.  Cheap: the waveform is already stored.
        """
        W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                      pss.irefnode, axis=0)
        if W.size == 0 or W.shape[1] < 4:
            return 0.0
        row = W[0]
        half = len(row) // 2
        den = float(np.max(np.abs(row)))
        if den <= 0.0:
            return 0.0
        return float(np.max(np.abs(row[:half] + row[half:2 * half]))) / den

    def _warn_if_orbit_is_asymmetric(self, pss):
        """⚠ On a strongly asymmetric orbit the modal sum carries an O(h)
        discretisation residual that the symmetric fixtures never show.

        The residual is the adjoint replay's `O(h)` discretisation error
        (~6 % at 400 points per period, halving per doubling, against a
        Monte-Carlo-validated Lyapunov reference; 1.0001 on a symmetric
        orbit), converging to 1 -- not a defect.  So this warns that the
        residual is grid-limited on such an orbit and says how to shrink it.

        History: `doc/shooting_history.md`,
        `PAC._warn_if_orbit_is_asymmetric`.
        """
        try:
            asym = self._orbit_asymmetry(pss)
        except Exception:
            return
        if asym <= self.ORBITAL_ASYMMETRY_LIMIT:
            return
        warnings.warn(
            'PAC.orbital_correlation: this orbit has half-wave asymmetry '
            '%.3f. On such an orbit the modal sum carries an O(h) '
            'discretisation residual of the adjoint replay that symmetric '
            'orbits do not show -- measured 6%% high at 400 points per period '
            'and halving per doubling against a Monte-Carlo-validated '
            'reference. Refine the grid to tighten it, or use '
            'PAC.oscillator_covariance (Lyapunov) for the covariance alone.'
            % (asym,),
            RuntimeWarning, stacklevel=3)

    def orbital_correlation(self, pss, H=None):
        """`R_yy(0)` and the `C_lhj` of Traversa & Bonani eq (22) — A9 step 3.

        Returns `(R, C)`.  `R` is the STATIONARY transverse (orbital)
        state covariance, `m x m` real symmetric — eq (23),
        `R = Σ_{l≥2,h,j} C_lhj`.  `C` maps `(l, h, j)` to the `m x m`
        complex coefficient, over every non-null orbital mode `l ≥ 2` and
        harmonics `|h|, |j|, |j'| ≤ H`.  `H` defaults to
        `ORBITAL_HARMONICS`; van der Pol converges by `H = 4`, and a
        strongly non-sinusoidal orbit needs more — check by raising it.

        ⚠ STATIONARY WHITE SOURCES ONLY, and that is what makes it
        computable without `B`.  Eq (22) needs the Fourier coefficients
        of `v_l(t)^T B(t)`; with `CY = B B^T` constant those products
        collapse to `V~_{l'k}^T CY V~*_{lk'}`, so the noise enters only
        through the reduced `CY` that `_cy_reduced` already refuses to
        hand over when it is bias-dependent.

        ⚠ `CY/2`, NOT `CY`.  The library's `CY` is one-sided; eq (22)
        integrates `B B^T` as a two-sided intensity -- the `kT/C`-calibrated
        Monte Carlo injection `Var(i) = CY/(2h)`, confirmed here by three
        routes agreeing.

        ⚠⚠ GATED THREE WAYS, because a modal sum transcribed from an image
        of an equation is exactly the object to distrust.  (i) This sum
        against `R_yy(0)` evaluated from its DEFINITION as a 1-D Lyapunov
        integral along the orbital mode, no Fourier machinery; (ii) both
        against the CYCLE-MEAN transverse part of `oscillator_covariance`'s
        samples, which shares no machinery with either.

        ⚠ THE REFERENCE IS THE CYCLE MEAN, NOT `K_orb(0)`.  Lemma 3.5's
        `R∞_yy` depends on `τ` only — the stationary part.  At `t = 0`
        van der Pol's amplitude direction is pure-v while this is
        isotropic, which is a rotating radial direction averaged over a
        cycle, not a disagreement.

        ⚠ AND THE REFERENCE IS OBLIQUELY PROJECTED.  Subtracting only the
        SECULAR growth `(t/T) d u u^T` from the Lyapunov samples leaves the
        phase direction's BOUNDED within-period variance, which eq (22)'s
        `l >= 2` sum correctly excludes.  Demir's `y` is defined by the
        OBLIQUE projection `v_1^T y = 0`: project the samples with
        `Pi = I - u v^T/(v^T u)`.  Unprojected, the reference is off by a
        residual falling as 1/Q.
        (The phase-orbital CORRELATION's tau = 0 value is not a missing
        term: eq (18a)'s brace is {1 - 1} = 0 there, and eq (23) states
        R_yy(0) = sum C_lhj alone.)

        History: `doc/shooting_history.md`, `PAC.orbital_correlation`.
        """
        self._refuse_coloured(pss, 'orbital_correlation')
        H = self.ORBITAL_HARMONICS if H is None else int(H)
        modes = pss.floquet_modes(pss)
        m = self.cir.n - 1
        Tp = float(pss.period)
        w0 = 2.0 * np.pi / Tp
        CY2 = 0.5 * np.real(np.asarray(self._cy_reduced(pss, 0.0)))
        ## the phase mode by its tangent alignment (see `_phase_mode_split`),
        ## and NEVER in the orbital sum: swept in with a near-zero exponent it
        ## blows up as 1/|mu|^2
        _kph, orb = self._phase_mode_split(pss, modes, 'PAC.orbital_correlation')
        if not orb:
            raise ValueError(
                'PAC.orbital_correlation: no orbital mode -- every non-null '
                'multiplier is the phase mode.')

        def fcoef(P):
            ## `_period_dft`: the index DFT on a uniform grid, unchanged, and
            ## the trapezoid-weighted sum at the true times on a non-uniform one
            X = np.asarray(P)[:, :-1]
            N = X.shape[1]
            return self._period_dft(pss, X.T).T, N

        U, V, N = {}, {}, None
        for k in orb:
            U[k], N = fcoef(modes[k]['p'])
            V[k], _ = fcoef(modes[k]['q'])
        H = min(H, N // 2 - 1)
        hs = np.arange(-H, H + 1)
        idx = lambda k: k % N

        C = {}
        R = np.zeros((m, m), dtype=complex)
        for l in orb:
            mul = modes[l]['mu']
            for lp in orb:
                mulp = modes[lp]['mu']
                for j in hs:
                    Ulj = U[l][:, idx(j)]
                    ## the Lambda products for every (h, j') at once
                    Vl_hj = V[l][:, idx(hs - j)]              # m x nh  (h - j)
                    for jp in hs:
                        res = 1.0 / (1j * (j - jp) * w0 - mulp - np.conj(mul))
                        outer = res * np.outer(U[lp][:, idx(jp)], np.conj(Ulj))
                        Vlp_hjp = V[lp][:, idx(hs - jp)]      # m x nh  (h - j')
                        sc = np.einsum('ih,ik,kh->h', Vlp_hjp, CY2, np.conj(Vl_hj))
                        for hi, h in enumerate(hs):
                            term = sc[hi] * outer
                            key = (l, int(h), int(j))
                            C[key] = C.get(key, 0.0) + term
                            R += term
        return np.real(R), C

    def orbital_spectrum(self, pss, offsets, output, harmonic=1, H=None):
        """`S_yy` — the ORBITAL (amplitude) noise spectrum. A9 step 4.

        Returns `S` at `harmonic*f0 + offsets`, in the same V^2/Hz scale as
        `oscillator_spectrum`'s `S_v`, so **the two are summed** — which is
        what Traversa & Bonani (TCAS-I 2011) say to do:

            x(t) = x_s(t + a(t)) + y(t)      a = phase, y = orbital

        with the phase--orbital CROSS term dropped.  ⚠ That is a documented
        approximation with a KNOWN SIGN, not an oversight: the paper reports
        the correlation spectrum negligible on two circuits, and that when
        present it *decreases* the total.  **Dropping it therefore OVER-states
        noise** — conservative for a design margin, wrong in a known
        direction.  It is identically zero with no AM-to-PM coupling.

        ⚠⚠ "NEGLIGIBLE" IS A PROPERTY OF THOSE CIRCUITS, and the
        over-statement can be large.  Against pnoise (the total linear
        sideband noise, confirmed by a Monte Carlo of the SDE), on van der
        Pol with an `a u^2` asymmetry, the sum is right to 0.1 % on a
        SYMMETRIC orbit and over-states by up to 3.2x (5 dB) at half-wave
        asymmetry 0.1, with no grid dependence.

        ⚠⚠ THE DROPPED CROSS TERM IS THE CAUSE -- WITH EVERY HARMONIC KEPT.
        `S_corr` from eq (92) keeps only the PPV's DC harmonic AT THE NOISE
        SOURCE's row (~1e-8 of the total where a tank inductor shorts that
        node at DC); the full-harmonic correlation is -1.1 to -2.4x this
        spectrum at `a = 0.30`, and `modal_spectrum`'s three terms sum to
        pnoise.  The over-statement is in the decomposition's
        frequency-independent terms above f_amp: the phase half is the
        Lorentzian's frequency-independent PPV, and the orbital half
        over-states by a factor FLAT in offset -- the orbital mode's AM
        share at the output, `sin^2 arg(U_{l,1}/U_{0,1})`.  Traversa &
        Bonani's own Figs 1-2 show the same limit.  A warning fires above
        `ORBITAL_ASYMMETRY_LIMIT`.

        **Lemma 3.5**: the orbital spectrum is a sum of Lorentzians centred at
        `j*w0 + Im(mu_l)` with half-width `|Re(mu_l)| + (1/2) h^2 w0^2 c`,
        weighted by the `C_lhj` of eq (22).  Every input already exists:
        `orbital_correlation` returns `C_lhj` (gated three ways), the
        exponents come from `floquet_modes`, and `c` from
        `diffusion_constant`.

        ⚠ UNIT CONVERSION, DONE ONCE HERE.  Lemma 3.5's widths are ANGULAR.
        `(1/2) h^2 w0^2 c` rad/s is `pi h^2 f0^2 c` Hz -- exactly the
        half-width `lorentzian` already uses for the phase line -- and
        `|Re(mu_l)|` rad/s is `|Re(mu_l)|/(2 pi)` Hz.  The two half-widths
        ADD, so an orbital mode's line is the phase line broadened by the
        mode's own relaxation rate.

        ⚠⚠ AND THAT IS WHY IT MATTERS AT LARGE OFFSET, WHICH IS THE WHOLE
        POINT OF THE ITEM.  The phase line's width is `pi h^2 f0^2 c`, which
        for a good oscillator is tiny, so its skirt has fallen as `1/f^2` long
        before the orbital line -- width `|Re(mu_2)|/(2 pi)`, i.e. the
        AMPLITUDE RELAXATION RATE -- has even started to roll off.  The
        crossover therefore sits near

            f_amp = -ln(lam2) f0 / (2 pi) = f0 / (2 pi Q)

        the same pole `oscillator_spectrum` warns about from the other side.
        ⚠ Those two arrived independently -- one from a commercial
        simulator's excess over our phase-only answer, one from this paper's
        modal sum -- and they must land in the same place.  That is the gate
        (`test_the_orbital_spectrum_crosses_the_phase_spectrum_near_f_amp`),
        and it is the check that can actually fail.

        ⚠ `output` follows `oscillator_spectrum`: an integer indexes the
        REDUCED state (the reference row already removed), an array is a
        weight vector over it.

        ⚠ STATIONARY WHITE SOURCES ONLY -- inherited from
        `orbital_correlation`, which needs `CY` constant for eq (22)'s
        products to collapse.

        History: `doc/shooting_history.md`, `PAC.orbital_spectrum`.
        """
        try:
            _asym = self._orbit_asymmetry(pss)
        except Exception:
            _asym = 0.0
        if _asym > self.ORBITAL_ASYMMETRY_LIMIT:
            ## ⚠ NOT the grid residual `_warn_if_orbit_is_asymmetric` names:
            ## the sum this spectrum is meant for over-states the TOTAL on an
            ## asymmetric orbit, and no refinement changes it -- see the
            ## docstring.
            warnings.warn(
                'PAC.orbital_spectrum: this orbit has half-wave asymmetry '
                '%.3f. On an asymmetric orbit S_ph + S_orb over-states the '
                'total sideband noise above f_amp: x1.03 / x1.45 / x3.2 at '
                'asymmetry 0.033 / 0.067 / 0.100 (van der Pol, C=4, Q=8, '
                '10 f_amp), confirmed by Monte Carlo, which agrees with '
                'pnoise to 1 %%. Refining the grid does not change it. The '
                'cause is the phase-orbital correlation this sum drops (with '
                'every harmonic kept it is -1.1 to -2.4x the orbital term '
                'there): use PAC.modal_spectrum for a phase/orbital/'
                'correlation split that sums to the total, or PAC.pnoise for '
                'the total.'
                % (_asym,),
                RuntimeWarning, stacklevel=2)
        R, C = self.orbital_correlation(pss, H=H)
        modes = pss.floquet_modes(pss)
        c = float(self.diffusion_constant(pss))
        f0 = 1.0 / float(pss.period)
        m = pss.cir.n - 1

        row = _output_row(output, m)

        ## ⚠⚠ NO ORBITAL LINE AT THIS HARMONIC -- refuse rather than return the
        ## tails of the others.  The line weight at `j f0` is
        ## `W_j = sum_{l,h} Re(row C_lhj row)`; where it is zero (a symmetric
        ## orbit's even harmonics: the modes' own Fourier content vanishes)
        ## what this would return is the neighbouring lines' Lorentzian tails
        ## (3.2x LOW at 2 f0 on van der Pol against a Monte Carlo).  ⚠ NOT
        ## caught: DC, where a small line can exist and the model reads ~100x
        ## HIGH (the tank inductor shorts the node, which Lorentzian tails do
        ## not know), and 2 f0 on an asymmetric orbit -- away from the
        ## fundamental use `pnoise`.
        _W = {}
        for (_l, _h, _j), _cl in C.items():
            _W[_j] = _W.get(_j, 0.0) + float(np.real(row @ _cl @ row))
        _Wtot = sum(abs(v_) for v_ in _W.values())
        if abs(_W.get(int(harmonic), 0.0)) <= 1e-9 * max(_Wtot, 1e-300):
            raise ValueError(
                'PAC.orbital_spectrum: no orbital line at harmonic %d for this '
                'output (weight %.3e of %.3e), so the value here would be the '
                'tails of other lines, not the noise at that frequency '
                '(measured 3.2x low at 2 f0 on a symmetric van der Pol). Use '
                'PAC.pnoise there.' % (int(harmonic), _W.get(int(harmonic), 0.0),
                                       _Wtot))

        f = float(harmonic) * f0 + np.atleast_1d(
            np.asarray(offsets, dtype=float))
        S = np.zeros_like(f, dtype=float)
        for (l, h, j), Clhj in C.items():
            ## The weight is the output's own share of this term.  It is real
            ## for the total (`R` is real symmetric); an individual `(l,h,j)`
            ## can carry a small imaginary part that cancels against its
            ## conjugate partner, so take the real part per term rather than
            ## asserting each is real.
            w = float(np.real(row @ Clhj @ row))
            if w == 0.0:
                continue
            mul = modes[l]['mu']
            ## Hz, both terms -- see the unit note above.
            gam = abs(float(np.real(mul))) / (2.0 * np.pi) \
                + np.pi * float(h) ** 2 * f0 ** 2 * c
            fc = float(j) * f0 + float(np.imag(mul)) / (2.0 * np.pi)
            if gam <= 0.0:
                continue
            ## Normalised Lorentzian: integrates to 1 over all `f`, so the
            ## total power is `sum(w) = row^T R row` by construction.
            S = S + w * (gam / np.pi) / ((f - fc) ** 2 + gam ** 2)
        return S

    def modal_spectrum(self, pss, offsets, output, harmonic=1, H=None,
                       sidebands=None):
        """Phase, orbital AND phase-orbital CORRELATION spectra from ONE modal
        transfer, which sum to the total.

        Returns a dict of arrays at `harmonic*f0 + offsets` (a negative offset
        is the lower sideband), on the scale of `oscillator_spectrum`'s `S_v`
        and `orbital_spectrum` (0.5x a one-sided PSD):

            'phase', 'orbital', 'correlation', 'total'
            total = phase + orbital + correlation

        Every Floquet mode `l` -- the phase mode (`mu = 0`) and each orbital
        mode -- carries noise from input sideband `m` to the output at `w`:

            T_m^l(w) = sum_j (d . U_{l,j}) V_{l,m-j}^T / (i(w - j w0) - mu_l + a_j)

        with `U`, `V` the Fourier coefficients of `p_l`, `q_l` (the
        conventions of `orbital_correlation`) and `a_j = j^2 w0^2 c / 2` the
        phase-diffusion rate of output harmonic `j`.  With `T = T^0 + sum_l
        T^l` the output is `sum_m T_m (CY/2) T_m^H`; `phase`, `orbital` and
        `correlation` are its phase-phase, orbital-orbital and
        `2 Re(phase-orbital)` blocks.

        ⚠⚠ WHY IT EXISTS.  `oscillator_spectrum(frequency_aware=False) +
        orbital_spectrum` over-states an asymmetric orbit's total by up to
        3.2x (Monte-Carlo-confirmed).  The missing piece IS the correlation
        -- but with every harmonic kept, not in the form Traversa & Bonani
        keep: their eq (92) retains only its DC harmonic, ~1e-8 of the total
        on van der Pol (the tank inductor shorts the source node at DC).  It
        removes the DC-PPV phase excess above f_amp AND the orbital mode's PM
        projection at the output -- the orbital line's AM share
        `sin^2 arg(U_{l,1}/U_{0,1})` is the flat factor `orbital_spectrum`
        over-states by.  `total / (pnoise/2)` is 1.0006 on a symmetric van
        der Pol; on an asymmetric one the excess is the modes' O(h) grid
        error (it halves with the grid).
        ⚠ Because the correlation cancels most of the other two, a few
        percent of error in any part is AMPLIFIED in the total -- which is
        why the three are computed together here rather than the
        correlation being offered as an add-on to `oscillator_spectrum +
        orbital_spectrum`: those line-shape spectra keep only the resonant
        term of each line (2.5 % short at 10 f_amp even on a symmetric
        orbit), which is harmless alone and not under cancellation.  The
        modal sum also reproduces pnoise's upper/lower sideband asymmetry,
        which the two-term sum cannot.

        Near the carrier `phase` IS the library Lorentzian (symmetric orbit)
        and `correlation` is ~1e-6 of it.  Above f_amp `total` agrees with
        `pnoise`; within the linewidth pnoise has no meaning and this is the
        route.

        ⚠ Stationary WHITE sources, free-running oscillators, and the dense
        `floquet_modes` only (inherited).  `harmonic >= 1`: harmonic 0 was
        never measured.  `H` defaults to `ORBITAL_HARMONICS` (capped by the
        grid), `sidebands` to `2 H`.  `output` follows `orbital_spectrum`.

        History: `doc/shooting_history.md`, `PAC.modal_spectrum`.
        """
        self._check_circuit(pss)
        self._refuse_coloured(pss, 'modal_spectrum')
        self._refuse_driven(pss, 'modal_spectrum')
        if int(harmonic) < 1:
            raise ValueError(
                'PAC.modal_spectrum: harmonic must be >= 1 -- harmonic 0 was '
                'never measured against pnoise. Use PAC.pnoise there.')
        modes = pss.floquet_modes(pss)
        ## the phase mode by its tangent alignment, on any grid -- see
        ## `_phase_mode_split` (a 1e-6 window on |lam| - 1 refuses every gear
        ## solve on a non-uniform grid)
        _kph, orb = self._phase_mode_split(pss, modes, 'PAC.modal_spectrum')
        ph = [_kph]
        m = pss.cir.n - 1
        row = _output_row(output, m)
        c = float(self.diffusion_constant(pss))
        w0 = 2.0 * np.pi / float(pss.period)
        CY2 = 0.5 * np.real(np.asarray(self._cy_reduced(pss, 0.0)))
        N = np.asarray(modes[ph[0]]['p']).shape[1] - 1
        H = self.ORBITAL_HARMONICS if H is None else int(H)
        H = min(H, N // 2 - 1)
        M = 2 * H if sidebands is None else int(sidebands)
        js = np.arange(-H, H + 1)
        ms = np.arange(-M, M + 1)
        a_j = 0.5 * js.astype(float) ** 2 * w0 ** 2 * c

        ## per mode: the output's share of each harmonic of p_l, and q_l's
        ## Fourier coefficients (m x N).  ⚠ The phase mode's exponent is set
        ## to 0 exactly: its multiplier is 1 to rounding, and a 1e-16 real part
        ## would put a spurious pole width on the Lorentzian.
        coef = []
        for l in ph + orb:
            ## ⚠ `_period_dft`, not an index DFT, which is 8-13 % off and does
            ## not converge on a 3:1 grid (see `PSS._period_quadrature`)
            Ul = self._period_dft(pss, np.asarray(modes[l]['p'])[:, :-1].T).T
            Vl = self._period_dft(pss, np.asarray(modes[l]['q'])[:, :-1].T).T
            coef.append((l, row @ Ul[:, js % N], Vl,
                         0.0 if l == ph[0] else complex(modes[l]['mu'])))

        def transfer(w, entry):
            ## (2M+1) x m; looped over j so memory stays m x (2M+1) per mode
            _l, u, Vl, mul = entry
            g = u / (1j * (w - js * w0) - mul + a_j)
            T = np.zeros((ms.size, m), dtype=complex)
            for ji, j in enumerate(js):
                if g[ji] != 0.0:
                    T += g[ji] * Vl[:, (ms - j) % N].T
            return T

        def quad(A, B):
            return complex(np.einsum('mi,ik,mk->', A, CY2, np.conj(B)))

        offs = np.atleast_1d(np.asarray(offsets, dtype=float))
        res = {k: np.zeros(offs.shape, dtype=float)
               for k in ('phase', 'orbital', 'correlation', 'total')}
        for i, o in enumerate(offs.ravel()):
            w = float(harmonic) * w0 + 2.0 * np.pi * float(o)
            Tp = transfer(w, coef[0])
            To = np.zeros_like(Tp)
            for entry in coef[1:]:
                To += transfer(w, entry)
            sp = float(np.real(quad(Tp, Tp)))
            so = float(np.real(quad(To, To)))
            sc = 2.0 * float(np.real(quad(Tp, To)))
            ix = np.unravel_index(i, offs.shape)
            res['phase'][ix] = sp
            res['orbital'][ix] = so
            res['correlation'][ix] = sc
            res['total'][ix] = sp + so + sc
        return res

    def correlation_spectrum(self, pss, offsets, output, harmonic=1, H=None,
                             sidebands=None):
        """The FULL-harmonic phase-orbital correlation spectrum:
        `modal_spectrum(...)['correlation']`.

        ⚠ It sums to the total with `modal_spectrum`'s own `phase` and
        `orbital` -- NOT with `oscillator_spectrum + orbital_spectrum`, whose
        line-shape approximations are a few percent off under the cancellation
        this term produces (see `modal_spectrum`).  Negative where AM-to-PM
        coupling exists: -1.1 to -2.4x the orbital term on van der Pol at
        half-wave asymmetry 0.10.
        """
        return self.modal_spectrum(pss, offsets, output, harmonic=harmonic,
                                   H=H, sidebands=sidebands)['correlation']

    #: the phase mode may sit this far off the unit circle before it is refused
    PHASE_MODE_MAX_DEPARTURE = 1e-3

    def _phase_mode_split(self, pss, modes, where):
        """`(phase_index, orbital_indices)` -- the phase mode identified by WHAT
        DEFINES IT, its eigenvector's alignment with the orbit tangent, not by
        a window on `|lam| - 1`.

        On a uniform grid the phase multiplier is 1 to rounding.  On a grid
        whose step varies, a multistep or trapezoidal solve loses time-
        translation symmetry and the multiplier leaves the circle at O(h^2)
        (radau keeps it at 1 to ~1e-11), so a window on `|lam| - 1` would
        refuse gear there.  The right eigenvector of the phase mode is the
        tangent `xdot(0)` (`C xdot = -i(x)` for the autonomous circuit),
        which no orbital mode shares, so alignment picks it on any grid;
        its exponent is then forced to 0 exactly, as the consumers already
        do.  The departure is WARNED with its size when it exceeds rounding,
        and the split is REFUSED when a second multiplier lies within ten
        times that departure of the circle with any alignment -- the case a
        window ever protected against.

        History: `doc/shooting_history.md`, `PAC._phase_mode_split`.
        """
        m = pss.cir.n - 1
        irn = pss.irefnode
        xr = np.asarray(pss.waveform[1], dtype=float)[:, 0]
        xf = np.concatenate((xr[:irn], np.zeros(1), xr[irn:]))
        i_red = np.delete(np.asarray(pss.cir.i(xf, pss.epar), dtype=float).ravel(), irn)
        C0 = np.asarray(pss._C_at(xr), dtype=float)
        try:
            xdot0 = np.linalg.solve(C0, -i_red)
        except np.linalg.LinAlgError:
            xdot0 = np.linalg.lstsq(C0, -i_red, rcond=None)[0]
        nx = float(np.linalg.norm(xdot0))
        dep, cos = [], []
        for md in modes:
            u = np.asarray(md['u0'])[:m]
            dep.append(abs(abs(complex(md['lam'])) - 1.0))
            cos.append(abs(complex(np.vdot(u, xdot0))) / max(float(np.linalg.norm(u)) * nx, 1e-300))
        cand = [k for k in range(len(modes)) if dep[k] <= self.PHASE_MODE_MAX_DEPARTURE and cos[k] > 0.9]
        if not cand:
            raise ValueError(
                '%s: no Floquet mode is both within %.0e of the unit circle and '
                'aligned with the orbit tangent (|lam|-1: %s; alignment: %s).  Is '
                'this an autonomous oscillator solved at its own period?'
                % (where, self.PHASE_MODE_MAX_DEPARTURE,
                   ', '.join('%.1e' % d for d in dep), ', '.join('%.2f' % c for c in cos)))
        k = max(cand, key=lambda j: cos[j])
        rival = [j for j in range(len(modes)) if j != k and dep[j] <= max(10.0 * dep[k], 1e-6)]
        if rival:
            raise ValueError(
                '%s: a second Floquet multiplier (|lam|-1 = %s) lies as close to '
                'the unit circle as the phase mode (%.1e), so the phase mode '
                'cannot be told from an orbital one.  A uniform grid or '
                'method=\'radau\' puts the phase multiplier at 1 to rounding.'
                % (where, ', '.join('%.1e' % dep[j] for j in rival), dep[k]))
        if dep[k] > 1e-6:
            warnings.warn(
                '%s: the phase mode sits %.1e off the unit circle (alignment '
                'with the orbit tangent %.4f); its exponent is forced to 0.  '
                'On a non-uniform grid a multistep or trapezoidal solve leaves '
                'it there at O(h^2); the modal parts are then the method\'s '
                'order (measured second order for gear on a 3:1 grid).'
                % (where, dep[k], cos[k]), RuntimeWarning, stacklevel=3)
        return k, [j for j in range(len(modes)) if j != k]

    def diffusion_constant(self, pss):
        """`c` — the phase diffusion constant, in seconds.

        `c = (1/T) ∫ v₁ᵀ(t) B(t) Bᵀ(t) v₁(t) dt` with `B Bᵀ = CY`, so this
        is the time-average of a QUADRATIC form in the PPV.  It is the one
        scalar the whole free-running phase-noise spectrum is built from,
        and it reads, for a designer, as JITTER PER SECOND.

        ⚠ QUADRATIC FOR WHITE SOURCES, LINEAR FOR COLOURED ONES, and the
        two are different functionals of the same vector: a coloured
        source contributes `V_0m = (1/T) ∫ v₁ᵀ B_cm dt`, with no square.
        Using this one for a coloured source returns a plausible non-zero
        number from the same PPV.  Only stationary white sources are
        supported here, which `_cy_reduced` enforces.

        ⚠ `CY/2`, as in `covariance`: settled against `kT/C`, which is
        external to both (an injection of `Var(i) = CY/h` per step
        reproduces 1.92x `kT/C`).  A Monte Carlo built on the convention
        under test cannot test it.

        ⚠ `ppv()` normalises on the FIRST BLOCK (`v[:m] . xdot = 1`), which
        is right for a perturbation entering the first block -- an injected
        current, and what every shipped path does -- but wrong for
        contracting against a full PAIR deviation, where the factor is
        `1/(v . u_pair)`.  The sign is a convention: a later zero crossing
        means DELAYED, while projecting onto the tangent makes positive
        mean ADVANCED.

        History: `doc/shooting_history.md`, `PAC.diffusion_constant`.
        """
        self._check_circuit(pss)
        self._refuse_coloured(pss, 'diffusion_constant')
        self._refuse_driven(pss, 'diffusion_constant')
        return self._white_diffusion_at(pss, 2.0 * np.pi / float(pss.period))

    def _refuse_driven(self, pss, what):
        if not getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.%s: phase diffusion is a property of a '
                "FREE-RUNNING oscillator. A driven circuit's phase is its "
                "source's, and its noise is pnoise's problem, not this one."
                % what)

    @staticmethod
    def _period_weights(tms, nsamp, T, pss=None):
        """Periodic TRAPEZOID weights for samples at `tms[0..nsamp-1]` over
        a period `T`: `w_j = (g_j + g_{j-1}) / 2` with `g_j` the gap to the
        next sample and the last gap closing the period.

        ⚠ NOT `h = diff(times)`: that is the LEFT RECTANGLE rule, spectrally
        accurate on a uniform periodic grid but FIRST order on a NON-UNIFORM
        one (its error is (1/2) integral h'(t) y(t) dt, not zero).  On a
        uniform grid `0.5 h + 0.5 h == h` exactly.  ⚠ AND THE TRAPEZOID IS
        ITSELF A SECOND-ORDER CAP on a smoothly varying grid: with `pss`
        given and a non-uniform grid, these are the periodic cubic-spline
        weights of `periodic_spline_weights`, broken at the nodes of landed
        events; a uniform grid is unchanged.

        History: `doc/shooting_history.md`, `PAC._period_weights`."""
        tms = np.asarray(tms, dtype=float).ravel()
        n = int(nsamp)
        t = tms[:n]
        g = np.empty(n)
        g[:-1] = t[1:] - t[:-1]
        g[-1] = float(T) + t[0] - t[-1]
        if (pss is not None and n >= 4
                and float(np.max(g)) / float(np.min(g)) - 1.0 > pss.UNIFORM_GRID_TOL):
            return periodic_spline_weights(t, T, pss._event_nodes(t, T))
        w = 0.5 * g
        w[1:] += 0.5 * g[:-1]
        w[0] += 0.5 * g[-1]
        return w

    def _white_diffusion_at(self, pss, w):
        """`(1/T) integral v_1^T (CY(w)/2) v_1 dt` with `CY` FROZEN at `w`.

        The white functional at one frequency, with no refusal: it is `c`
        when the source is white, and for a coloured source it is the
        value `diffusion_constant` refuses.  `phase_psd`
        reads it at the carrier for the Lorentzian CORNER, which is a
        white-noise construct whatever the source's colour; the spectrum
        itself comes from `coloured_diffusion_resolved`.

        History: `doc/shooting_history.md`, `PAC._white_diffusion_at`.
        """
        v, info = pss.ppv()
        ## `lambda_2` is computed here anyway; `oscillator_spectrum` needs it to
        ## report its own validity limit and a second `ppv()` would be a full
        ## extra solve.  Recorded, not returned, so this method's signature is
        ## unchanged -- and read ONLY immediately after a call, which is how
        ## `oscillator_spectrum` uses it.
        self._last_second_multiplier = (
            info.get('second_multiplier'),
            info.get('second_multiplier_certified'))
        m = pss.cir.n - 1
        ## ⚠ `samples_eq`, NOT `samples`.  `CY` is an EQUATION-ROW
        ## covariance and `samples` is `C^T v_1`; contracting that would make
        ## `c` wrong by `C^2` on the differential rows and exactly zero on
        ## the algebraic ones.  See `_equation_row_ppv`.
        S = np.asarray(info['samples_eq'])[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        ## the samples' own orbit, not `pss.period` -- see `ppv()`'s 'period'
        T = float(info['period'])
        h = self._period_weights(tms, S.shape[0], T, pss)
        cy = self._cy_reduced(pss, float(w))
        ## ⚠ A NOISE SOURCE ON AN INDEX-2 CONSTRAINT GIVES c = 0, SILENTLY: a
        ## voltage noise in series with a DC source inside a capacitor loop
        ## perturbs an algebraic constraint -- a DIFFERENTIATED input, whose
        ## response is a charge jump the PPV projection cannot represent --
        ## and its share of `c` is exactly 0 for every method.  Named here
        ## once: the PPV's algebraic fallback fired (index >= 2) and `CY` has
        ## power on an algebraic row.  (The gate is the index-2 condition
        ## itself -- `G[A,Z]` singular at the orbit point, the test the PPV's
        ## algebraic fallback makes -- computed here so every kind is
        ## covered.)
        try:
            x0r = np.asarray(pss._period_state[1], dtype=float).ravel()
            irn = pss.irefnode
            x0f = np.concatenate((x0r[:irn], np.zeros(1), x0r[irn:]))
            _arows, _acols = pss._algebraic_adjoint_pattern(x0f)
            _idx2 = False
            if _arows and _acols and len(_arows) == len(_acols):
                _Gz = np.asarray(pss._G_at(x0r), dtype=float)[np.ix_(
                    np.asarray(_arows, dtype=int), np.asarray(_acols, dtype=int))]
                _sv = np.linalg.svd(_Gz, compute_uv=False)
                _idx2 = (float(_sv[-1]) <= 1e-12 * max(float(_sv[0]), 1e-300))
            elif _arows:
                _idx2 = True
        except Exception:                                      # noqa: BLE001
            _idx2 = False
        if _idx2:
            try:
                _dcy = np.abs(np.real(np.diag(cy)))
                _on_alg = [r for r in _arows if _dcy[r] > 0.0]
                if _on_alg and float(np.max(_dcy)) > 0.0:
                    warnings.warn(
                        'PAC.diffusion_constant: this circuit is index >= 2 and '
                        'a noise source sits on an algebraic row (%s). A '
                        'perturbation of an index-2 constraint is a '
                        'differentiated input -- a charge jump -- which the '
                        'PPV projection cannot represent, and its share of c '
                        'is 0 here whatever the source (measured: exactly 0 '
                        'for every method). Only the differential rows\' '
                        'sources are counted.' % (_on_alg,),
                        RuntimeWarning, stacklevel=3)
            except Exception:                                  # noqa: BLE001
                pass
        ## ⚠ `cy/2`, THE SAME ONE-SIDED-TO-TWO-SIDED CONVERSION `covariance`
        ## USES.  `CY` is a one-sided density (a resistor's `4kT/R`); an
        ## injection of `Var(i) = CY/h` per step reproduces `1.92x kT/C`, so
        ## that convention carries TWICE the physical noise power.
        quad = np.einsum('ij,jk,ik->i', S, 0.5 * np.real(cy), S)
        return float((quad * h).sum() / T)

    def colour_projection(self, pss):
        """`<v_1>` — the PPV's TIME AVERAGE, which is a different functional.

        Returns `(vbar, info)`.  `vbar` is `(1/T) integral v_1(t) dt` over
        the orbit; `info` carries the per-entry `rms` and the ratio
        `|mean|/rms`, which is the number that says whether a coloured
        source at that node can upconvert at all.

        ⚠ COLOURED SOURCES CONTRACT THE SQUARE OF THE MEAN; WHITE ONES
        CONTRACT THE MEAN OF THE SQUARE.  `diffusion_constant` computes
        `(1/T) integral v^T (CY/2) v dt`.  A coloured source's low-frequency
        power cannot be modulated away, so what survives is
        `V_0m = (1/T) integral v_1^T B_cm dt` -- LINEAR, no square -- and the
        contraction is `vbar^T (CY/2) vbar`.  Same vector, same matrix,
        the mean and the square exchanged.

        ⚠ AND USING THE QUADRATIC ONE FOR A COLOURED SOURCE RETURNS A
        PLAUSIBLE NUMBER, NOT AN ERROR.  It is never zero where the white
        answer is not, so nothing downstream would look wrong -- and the
        two are not close approximations of each other (22 orders apart on
        van der Pol), so they cannot be substituted.

        ⚠ TWO INDEPENDENT MECHANISMS FORCE `vbar` TO ZERO, AND ONLY ONE OF
        THEM IS THE ONE DESIGNERS KNOW.  NEITHER ASYMMETRY ALONE NOR LOSS
        ALONE UPCONVERTS (on an LC oscillator with an even nonlinearity
        term and a series tank resistance, `Gamma/c` ~ 1e-22 with either
        alone and 1e-4 .. 1e-3 with both, while `c` barely moves):

        * A SYMMETRIC waveform gives `vbar = 0` -- Hajimiri & Lee, and the
          reason symmetry is the first thing a VCO designer reaches for.
        * A LOSSLESS LC TANK gives `vbar[0] = 0` STRUCTURALLY, whatever the
          waveform does.  `v` behaves as `C^T v_1` and `dv/dt = G^T v_1`,
          whose inductor row is exactly `v[0]`; periodicity of `v[1]` then
          forces `integral v[0] dt = 0`.  ⚠ THIS IS A PROPERTY OF THE
          TOPOLOGY, NOT OF THE ORBIT, and it is why van der Pol reports
          zero at every asymmetry -- it makes van der Pol useless as a
          POSITIVE fixture and perfect as a negative one.

        ⚠ `Gamma <= c` ALWAYS, at the same `CY`, by Cauchy-Schwarz on the
        weighted mean -- with equality only if `v` is constant over the
        orbit.  Both use the same quadrature here so the bound holds
        exactly at the discrete level, which makes it an assertion rather
        than an expectation.

        History: `doc/shooting_history.md`, `PAC.colour_projection`.
        """
        self._check_circuit(pss)
        if not getattr(pss, 'autonomous', False):
            raise ValueError(
                'PAC.colour_projection: the PPV time-average is the kernel '
                "of a FREE-RUNNING oscillator's coloured-noise upconversion. "
                "A driven circuit's phase is its source's.")
        _v, info = pss.ppv()
        m = pss.cir.n - 1
        ## ⚠ the EQUATION-ROW adjoint, for the same reason
        ## `diffusion_constant` uses it: a coloured source is an
        ## equation-row input too.
        S = np.asarray(info['samples_eq'])[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        T = float(info['period'])
        h = self._period_weights(tms, S.shape[0], T, pss)
        ## ⚠ THE SAME QUADRATURE `diffusion_constant` USES, deliberately:
        ## it is what makes `Gamma <= c` exact rather than approximate.
        vbar = (S * h[:, None]).sum(0) / T
        rms = np.sqrt((S ** 2 * h[:, None]).sum(0) / T)
        with np.errstate(divide='ignore', invalid='ignore'):
            sym = np.where(rms > 0, np.abs(vbar) / rms, 0.0)
        return vbar, {'rms': rms, 'symmetry': sym,
                      'samples': S, 'times': tms}

    def coloured_diffusion(self, pss, freqs):
        """`Gamma(f) = vbar^T (CY(2 pi f)/2) vbar` — the coloured analogue of `c`.

        Returns an array over `freqs`.  `CY` is evaluated at each offset,
        so a source whose density varies with frequency -- which is what
        "coloured" means -- is folded in exactly as a white one is.

        ⚠ NO FILTER, NO EXTRA STATE, NO SDE.  Demir 1996 synthesises 1/f
        from white sources through a Lorentzian network at "one state
        variable per decade", because Ito theory admits only white driving
        noise.  That is an artefact of the SDE formulation.  This path
        never forms an SDE, so a coloured source is just a different
        `S(f)` -- a SLOPE, NOT A STATE.  A commercial RF simulator confirms by omission:
        no filter and no augmentation in its treatment of flicker.

        ⚠ THE `CY/2` IS THE SAME ONE-SIDED-TO-TWO-SIDED CONVERSION THE
        REST OF THIS CLASS USES, and it is shared rather than repeated so
        the functions cannot drift apart over that factor.

        History: `doc/shooting_history.md`, `PAC.coloured_diffusion`.
        """
        vbar, _ = self.colour_projection(pss)
        out = []
        for f in np.atleast_1d(np.asarray(freqs, dtype=float)):
            cy = np.real(self._cy_reduced(pss, 2.0 * np.pi * float(f)))
            out.append(float(vbar @ (0.5 * cy) @ vbar))
        return np.asarray(out)

    def coloured_diffusion_resolved(self, pss, freqs, harmonics=None):
        """`c(f) = sum_l V_l^H (CY(2 pi |f - l f_0|)/2) V_l` — the fold PER HARMONIC.

        `V_l` are the Fourier coefficients of the equation-row PPV `v_1(t)`
        (the rows `diffusion_constant` contracts), so a source's density is
        read at the SOURCE-SIDE frequency `f - l f_0` for each harmonic it
        folds through -- which is what `pnoise` has done from the start and
        what a coloured source requires.  Returns an array over `freqs`.

        ⚠ `c + Gamma(f)` IS NOT THIS: `c` reads `CY` at ONE frequency
        (`2 pi / T`) as if it held at every harmonic, and `Gamma` is exactly
        the `l = 0` term of this sum, so `c + Gamma` counts `l = 0` twice.
        Neither shows on van der Pol, whose PPV at the tank node averages to
        zero (the inductor shorts the node at DC, so no core can bias it).

        EXACT FOR WHITE, BY PARSEVAL: with `CY` constant the sum is
        `(1/T) integral v_1^T (CY/2) v_1 dt = c`, and the discrete version
        with the grid's step weights reproduces `diffusion_constant` to
        round-off -- that equality pins the transform's normalisation, and
        it is asserted.  For a DC-centred colour (Lorentzian, flicker) and
        `f << f_0` the `l != 0` terms read `CY(l f_0)` to `O(f/f_0)`, so
        the sum differs from `c + Gamma` only where `V_0` is not small.

        `harmonics` caps `|l|`; by default every harmonic carrying more
        than 1e-14 of the PPV's energy is kept, which is all of them that
        can move the sum at double precision.

        History: `doc/shooting_history.md`, `PAC.coloured_diffusion_resolved`.
        """
        self._check_circuit(pss)
        self._refuse_driven(pss, 'coloured_diffusion_resolved')
        m = pss.cir.n - 1
        v0, info = pss.ppv()
        S = np.asarray(info['samples_eq'], dtype=float)[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        n = S.shape[0]
        ## the samples' orbit: both the 1/T AND the harmonic frequency `w0`
        ## must come from it, or Parseval leaks (1.5e-09 under a trap twin)
        T = float(info['period'])
        ## ⚠ THE SAME QUADRATURE `diffusion_constant` USES: one sample per
        ## step, weighted by that step, so that Parseval closes exactly.
        t = tms[1:1 + n]
        h = self._period_weights(t, n, T, pss)
        w0 = 2.0 * np.pi / T
        L = n // 2 if harmonics is None else int(harmonics)
        ls = np.arange(-L, L + 1) if harmonics is not None else np.arange(-L, L)
        E = np.exp(-1j * np.outer(ls, w0 * t)) * h[None, :]          # (nl, n)
        V = (E @ S) / T                                               # (nl, m)
        energy = np.sum(np.abs(V) ** 2, axis=1)
        keep = energy > 1e-14 * energy.sum()
        ls, V = ls[keep], V[keep]
        out = []
        for f in np.atleast_1d(np.asarray(freqs, dtype=float)):
            tot = 0.0
            for l, vl in zip(ls, V):
                cy = np.real(self._cy_reduced(pss, 2.0 * np.pi * abs(float(f) - l / T)))
                tot += float(np.real(np.conj(vl) @ (0.5 * cy) @ vl))
            out.append(tot)
        return np.asarray(out)

    def phase_psd(self, pss, offsets, harmonic=1):
        """`S_phi(f)` in rad^2/Hz at `offsets` from harmonic `i` — white AND coloured.

            S_phi,i(f) = i^2 f_0^2 c(f) / f^2,   c(f) = sum_l V_l^H (CY(f - l f_0)/2) V_l

        `c(f)` is `coloured_diffusion_resolved`: the phase diffusion with
        each harmonic's colour read at its own source-side frequency.  For
        each harmonic's colour read at its own source-side frequency.  For
        a white source it is `c` exactly (`c + Gamma(f)` would count the
        `l = 0` term twice).

        ⚠ THE CONVENTION IS PINNED BY `oscillator_spectrum`, NOT ARGUED.
        `lorentzian`'s far skirt is `i^2 f_0^2 c / f^2` exactly, and that
        object was gated by power conservation to 1.000000.  So this
        expression is the same quantity its tail already reports, with the
        coloured term added -- no second convention is introduced.

        ⚠ THE TWO TERMS ADD BECAUSE THE SOURCES ARE INDEPENDENT, and with
        `CY ~ 1/f` the coloured term gives `S_phi ~ 1/f^3` -- Kundert's
        "S_u(f) is generally pink ... then S_phi(f) would be proportional
        to 1/f^3 at low frequencies".

        ⚠ AND THIS IS THE LINEARISED PHASE MODEL, WHICH IS EXACT ENOUGH
        ONLY BECAUSE THE SOURCES ARE STATIONARY.  Vanassche, Gielen &
        Sansen (ICCAD 2002) locate the split between the exact phase
        equation `theta' = eps Gamma(t + theta) n(t)` and the approximate
        `theta' = eps Gamma(t) n(t)`: for a STATIONARY source the two
        "will, up to 0-th order in eps, predict the same output phase
        noise", and they diverge otherwise.  Their operational form is
        better than "non-stationary" -- "at first, near t = 0, the
        predicted phases are the same. However, when THETA BECOMES TOO
        LARGE [they diverge]".  A stationary source makes `theta` DIFFUSE;
        a driven one makes it grow SECULARLY, which is what carries it out
        of range.  ⚠ So this is sound for free-running noise and must NOT
        be reused for injection locking, a PLL in lock, or coupled
        oscillators -- there the shift has to stay inside the argument.

        ⚠ REFUSED BELOW THE LORENTZIAN CORNER, and this is a validity
        boundary rather than a conditioning one.  There the excess phase is
        a Wiener process whose spectrum is singular at the origin; the
        finite value the real lineshape attains comes from the NONLINEAR
        phase-to-voltage map, which `oscillator_spectrum` carries and this
        does not.  Reporting `S_phi` near the carrier is the mistake this
        object invites, so it raises instead.

        History: `doc/shooting_history.md`, `PAC.phase_psd`.
        """
        self._check_circuit(pss)
        f0 = 1.0 / float(pss.period)
        i = int(harmonic)
        if i < 1:
            raise ValueError('PAC.phase_psd: harmonic must be >= 1.')
        offs = np.atleast_1d(np.asarray(offsets, dtype=float))
        if np.any(offs <= 0.0):
            raise ValueError(
                'PAC.phase_psd: offsets must be positive; S_phi diverges '
                'at zero offset and that divergence is physical.')
        cres = self.coloured_diffusion_resolved(pss, offs)
        ## ⚠ THE CORNER IS THE WHITE LORENTZIAN'S, read at the carrier.  For
        ## a coloured source `f_h = pi i^2 f0^2 c` is not a lineshape
        ## parameter at all -- there is no Lorentzian -- and taking the
        ## folded value nearest the carrier instead would put a 1/f source's
        ## corner ABOVE the offsets, in front of the power bound below, which
        ## is the floor that actually binds for colour.
        c = self._white_diffusion_at(pss, 2.0 * np.pi * f0)
        ## The i-th harmonic's Lorentzian half-width.  `S_i(f) =
        ## i^2 f0^2 c / (pi^2 i^4 f0^4 c^2 + f^2)` is a Lorentzian in `f`
        ## whose denominator is `f_h^2 + f^2`, so `f_h = pi i^2 f0^2 c`.
        corner = np.pi * (i ** 2) * (f0 ** 2) * c
        if offs.min() <= corner:
            raise ValueError(
                'PAC.phase_psd: offset %.6g Hz is at or below the '
                'Lorentzian corner %.6g Hz for harmonic %d, where S_phi is '
                'not the right object -- the excess phase is a Wiener '
                'process and its spectrum is singular at the origin. The '
                'finite value the LINESHAPE attains there comes from the '
                'nonlinear phase-to-voltage map: use oscillator_spectrum().'
                % (float(offs.min()), corner, i))
        sphi = (i ** 2) * (f0 ** 2) * cres / offs ** 2

        ## ⚠ POWER CONSERVATION AS A SECOND, INDEPENDENT FLOOR -- and for a
        ## COLOURED source it is the binding one, by orders.  The
        ## normalised lineshape integrates to 1, and the integral over one
        ## box of width `df` on each side is a lower bound on it, so
        ##
        ##     2 df S_phi(df) <= 1
        ##
        ## is NECESSARY for the linearised skirt to be consistent with
        ## unit power.  Vanassche, Gielen & Sansen (2003) derive the same
        ## statement for a 1/f input as `df_c >= eps f0 sqrt(2 f_1f)`; the
        ## form here needs no assumption about the source's colour and
        ## reproduces their worked example exactly.
        ##
        ## ⚠ THE LORENTZIAN CORNER ABOVE DOES NOT CATCH THIS.  It is built
        ## from `c` alone, so it knows nothing about a `Gamma(f)` that
        ## grows as the offset falls (on this class's own flicker fixture
        ## the power bound bites 306x above the Lorentzian corner).
        ##
        ## ⚠ IT IS A LOWER BOUND ON THE BREAKDOWN, NOT THE BREAKDOWN.
        ## Passing it is not a guarantee (on Vanassche's own example the
        ## observed flattening sits at 3x the bound): this refuses what is
        ## definitely invalid and admits a band that is already suspect --
        ## deliberately, because refusing at 3x would be fitting a
        ## threshold to one example.
        ## ⚠ AND THE DERIVATION HAS A PRECONDITION THE BOUND DOES NOT
        ## STATE, so it is checked rather than assumed.  The box argument
        ## is `2 df S(df) <= integral_{-df}^{+df} S <= 1`, and the FIRST
        ## inequality needs `S(f) >= S(df)` for every `|f| <= df` -- the
        ## spectrum must not dip below its edge value anywhere further in.
        ## True of a monotone skirt, of the flattened near-carrier shape,
        ## and even with a spur, which ADDS power inside.
        ##
        ## ⚠ FALSE FOR A LOCKED PLL, whose phase-noise transfer function
        ## is HIGH-PASS: the spectrum dips below its edge value everywhere
        ## inside the loop bandwidth.  The bound is then NO LONGER DERIVED,
        ## and a floor that is not derived cannot be used as one.  This
        ## method refuses a driven circuit today; the check is there for the
        ## driven-oscillator work.
        probe = np.unique(np.concatenate((
            offs, np.logspace(np.log10(offs.min() / 1e3),
                              np.log10(offs.max()), 32))))
        sprobe = ((i ** 2) * (f0 ** 2)
                  * self.coloured_diffusion_resolved(pss, probe) / probe ** 2)
        if np.any(np.diff(sprobe) > 1e-12 * np.abs(sprobe[:-1])):
            k = int(np.argmax(np.diff(sprobe) > 0)) + 1
            raise ValueError(
                'PAC.phase_psd: the spectrum RISES with offset near '
                '%.6g Hz, so it dips below its edge value further in and '
                'the power bound below is no longer derived -- its box '
                'argument needs S(f) >= S(df) for every |f| <= df. That '
                'happens for a high-pass-shaped spectrum such as a locked '
                'loop, and for a source whose density grows faster than '
                'f^2. The bound may still hold; it is not established '
                'here, so it is refused rather than applied.'
                % float(probe[k]))

        power = 2.0 * offs * sphi
        bad = power >= 1.0
        if np.any(bad):
            k = int(np.argmax(bad))
            raise ValueError(
                'PAC.phase_psd: at offset %.6g Hz the linearised skirt '
                'already carries %.3f times the TOTAL power of the '
                'carrier (2 f S_phi >= 1), so it has broken down there -- '
                'a normalised spectrum integrates to 1. This bound is '
                'independent of the Lorentzian corner (%.6g Hz here) and '
                'for a coloured source it binds far earlier, because '
                'Gamma(f) grows as the offset falls. Sweep above it, or '
                'use oscillator_spectrum() for the lineshape. Note the '
                'TRUE breakdown is higher still: this is a lower bound.'
                % (float(offs[k]), float(power[k]), corner))
        return sphi

    @staticmethod
    def lorentzian(offsets, c, f0, harmonic=1):
        """The `i`-th harmonic's normalised lineshape at `offsets` from it.

            S_i(f) = i² f₀² c / (π² i⁴ f₀⁴ c² + f²)

        ⚠ EXACT FOR WHITE SOURCES, not a limiting form.  With coloured
        sources the transform "does not have a simple closed form" and only
        two-regime approximations exist — which is why the diffusion
        constants that feed it (`diffusion_constant`,
        `frequency_aware_diffusion`) refuse a coloured source
        (`_refuse_coloured`).

        ⚠ AND ITS TOTAL POWER IS EXACTLY 1.  `∫ a/(b²+f²) df = aπ/b`, and
        here `a = i² f₀² c`, `b = π i⁴ f₀⁴ c² ^ ½`… concretely `b = π i²
        f₀² c`, so the integral is exactly one.  **The carrier's power is
        redistributed, never created or destroyed** — which is the
        invariant that separates this from LTV small-signal treatments,
        which "erroneously predict infinite noise power [at the carrier] as
        well as infinite total integrated power".  It is asserted in the
        suite.

        The half-width is `π i² f₀² c` and the peak `1/(π² i² f₀² c)`, so a
        higher harmonic has a skirt scaling as `i²` and a corner as `i⁴` —
        `20 log₁₀(i)` dB noisier far out.
        """
        i = int(harmonic)
        if i == 0:
            return np.zeros_like(np.asarray(offsets, dtype=float))
        f = np.asarray(offsets, dtype=float)
        a = (i * i) * f0 * f0 * c
        b = np.pi * (i * i) * f0 * f0 * c
        return a / (b * b + f * f)

    def band_spread(self, pss, output, band, points=9, harmonic=1,
                    quantity='pnoise', **kw):
        """How much `S(r)·r²` VARIES across a band — the number that says
        whether a band mean and a point value are the same measurement.

        Returns `(spread, info)` with `spread = max/min` of `S(r)·r²` over
        `points` offsets spanning `band = (r_lo, r_hi)` in units of `f0`,
        and `info` carrying the samples, the band MEAN, the value at the
        band's midpoint, and their ratio.

        ⚠⚠ WHY THIS EXISTS.  Far above the AM corner both AM and PM fall as
        `1/r²`, so `S·r²` is flat and a band mean IS a point value, which
        makes the distinction invisible.  It is NOT general: a source behind
        a slow RC node has an in-band spectrum that is not `1/r²` at all
        (its `k = 0` term is filtered at the RC corner while the `k >= 1`
        terms are not, and their mix moves across the band), so a band mean
        and a point value can differ by ~4 % -- larger than most of the
        agreements this file asserts.  A comparison that takes a band mean
        on one side and a point value on the other is then measuring the
        convention, not the physics.

        ⚠ So: call this before comparing a measured band-averaged number
        against a computed point value, or vice versa.  A spread near 1
        licenses the shortcut; anything else says put both sides on the same
        footing.  `quantity` selects the surface (`'pnoise'`, `'S_pm'`,
        `'S_am'`, `'oscillator_spectrum'`); `**kw` is forwarded to it.

        History: `doc/shooting_history.md`, `PAC.band_spread`.
        """
        import numpy as _np
        f0 = 1.0 / float(pss.period)
        rs = _np.linspace(float(band[0]), float(band[1]), int(points))
        vals = []
        for r in rs:
            f = float(r) * f0
            if quantity == 'oscillator_spectrum':
                Sv, _i = self.oscillator_spectrum(pss, _np.array([f]), output,
                                                  harmonic=harmonic)
                v = float(_np.real(Sv[0]))
            elif quantity in ('S_pm', 'S_am'):
                am, pm, _b = self.am_pm_noise(pss, f, output, carrier=harmonic,
                                              **kw)
                v = float(_np.real(pm if quantity == 'S_pm' else am))
            else:
                v = float(_np.real(self.pnoise(pss, f, output, **kw)[0]))
            vals.append(v * float(r) ** 2)
        vals = _np.asarray(vals, dtype=float)
        lo = float(_np.min(_np.abs(vals)))
        spread = float(_np.max(_np.abs(vals)) / lo) if lo > 0.0 else _np.inf
        mean = float(_np.mean(vals))
        mid = float(_np.interp(0.5 * (rs[0] + rs[-1]), rs, vals))
        return spread, {'offsets': rs, 'values': vals, 'band_mean': mean,
                        'midpoint': mid,
                        'mean_over_point': (mean / mid) if mid != 0.0 else _np.inf}

    def oscillator_spectrum(self, pss, offsets, output, harmonic=1,
                            frequency_aware=True):
        """Free-running output spectrum at `offsets` from harmonic `harmonic`.

        ⚠⚠ THIS DOES NOT GO THROUGH `pnoise`'s SIDEBAND FOLD, AND IT CANNOT.
        The fold is a FREQUENCY-CONVERSION computation, complete for a
        driven circuit and structurally incomplete for an AUTONOMOUS
        oscillator: what it omits is exactly the near-carrier phase-noise
        skirt this method returns.  Rizzoli, Mastri & Masotti (IEEE MTT
        42-807, 1994): "frequency-conversion techniques alone are not
        sufficient ... for general autonomous circuits", because the
        noise-induced FREQUENCY MODULATION OF THE CARRIER at low offsets is
        not a frequency-conversion effect.  Their Section III names the two
        stacks -- CONVERSION noise, rising as 1/f for f -> 0, and MODULATION
        noise, "a jitter of the oscillatory steady state", rising as 1/f^3
        -- which DECOUPLE exactly at the steady state, and are "usually
        nearly equal" in an intermediate offset band (a cross-stack
        agreement test not built here).  Diagnostic value: a FLAT PSD near
        the carrier is neither slope -- it is the Phi(T) - I singularity.

        So the Floquet/PPV stack (`ppv`, `diffusion_constant`, this method)
        and the sideband fold (`pnoise`) ARE NOT TWO IMPLEMENTATIONS OF ONE
        QUANTITY, and unifying them is not a simplification waiting to be
        made.  ⚠ THE HAZARD IS THAT THE WRONG ONE STILL RETURNS A NUMBER:
        oscillator phase noise from the fold alone is a spectrum missing the
        dominant contribution near the carrier.  That is the completeness
        argument for the split; the efficiency argument (Floquet is cheaper)
        is the weaker one.

        Returns `(S_v, L_dBc)`.  ⚠ `S_v` is the Lorentzian lineshape scaled by
        `|X_1|^2 = A^2/4`, the carrier PHASOR's square -- which is HALF the
        carrier power `A^2/2` a one-sided PSD carries, so `S_v` is exactly
        0.5000x a one-sided PSD of the output voltage (against a reference
        simulator at every offset over four decades).  The scale is kept
        rather than doubled because callers may already divide by `|X_1|^2`
        themselves.  `L_dBc` is `S_v` normalised to the harmonic's own
        power, in dBc/Hz, and is unaffected.

        ⚠ NO SWEEP AND NO PER-FREQUENCY SOLVE.  Once the PSS waveform's
        Fourier coefficients and the scalar `c` are known, "we have an
        analytical expression that gives us the spectrum at any frequency.
        The computation of the spectrum is not performed separately for
        every frequency of interest."  Which also means it never meets the
        near-carrier singularity that a swept small-signal computation
        would, and never meets the 1/f sweep-grid trap — there is no sweep
        to place a point on.

        ⚠⚠ SCOPE: A SOURCE BEHIND A SLOW NODE.  The DC-PPV Lorentzian, for
        a noise source that reaches the core through a slow path (RC leg,
        tau >> T), holds only BELOW the source's corner `T/(2 pi tau)`;
        above it the true skirt is scaled by the PPV-harmonic-weighted
        filter `sum_k |G_k|^2 F_k(f) / sum_k |G_k|^2 F_k(0)` (G_k the PPV
        entry's Fourier coefficients at the source node, F_k the path's
        transfer at k f0 + f) -- 1/1000 at 0.1 f0 on a one-RC-leg fixture.
        `c` is still right, and a Monte Carlo of `c` cannot see it.
        `frequency_aware=True` (the default) replaces `c` by `c(f)` from the
        frequency-aware PPV (`frequency_aware_diffusion`), which matches
        pnoise through the corner (test ..._behind_a_slow_node_...) and to
        <= 2 % above f_amp on an orbit with AM-to-PM coupling.
        `frequency_aware=False` is the closed form, one `c` for every offset,
        and costs no solve.

        ⚠ AND IT IS THE ONLY ROUTE THAT IS VALID BELOW THE CORNER.  A
        small-signal analysis cannot produce `L(f)` there however well
        conditioned it is: the excess phase is a Wiener process, its
        spectrum has a singularity at the origin and no physical meaning,
        and the finite value `L` attains comes from the NONLINEAR
        phase-to-voltage map — which is what this closed form carries.
        Reporting `S_phi` near the carrier instead is the mistake that
        object invites.

        History: `doc/shooting_history.md`, `PAC.oscillator_spectrum`.
        """
        c = self.diffusion_constant(pss)
        f0 = 1.0 / float(pss.period)
        self._warn_above_amplitude_pole(offsets, f0)
        X = self.carrier_phasor(pss, output, harmonic)
        ## ⚠⚠ NO CARRIER, NO LINE -- AND THE ANSWER WOULD BE A PLAUSIBLE ZERO.
        ## This is a LINE-SHAPE model: it broadens the carrier's own harmonic.
        ## Where the output has no component at `harmonic` (a half-wave
        ## symmetric orbit's even harmonics, or `harmonic = 0`, which the
        ## Lorentzian returns as zeros by construction) there is nothing to
        ## broaden, and the true density is BROADBAND noise this method does
        ## not represent (~0 against 5.6e-6 V^2/Hz at 2 f0 on a symmetric
        ## van der Pol).  Refused, as `am_pm` refuses the same case.  ⚠ NOT
        ## caught: away from the fundamental the model also misses where a
        ## line DOES exist (asymmetric orbit, 2 f0: 0.40 of a Monte Carlo)
        ## -- use `pnoise` away from the fundamental.
        _scale = float(np.max(np.abs(self._output_waveform_row(pss, output))))
        if int(harmonic) == 0 or abs(X) <= 1e-9 * max(_scale, 1e-300):
            raise ValueError(
                'PAC.oscillator_spectrum: the output carries no component at '
                'harmonic %d (|X| = %.3e against a signal scale of %.3e), so '
                'there is no line for this line-shape model to broaden; the '
                'noise there is broadband and this method would return ~0 '
                '(measured: ~0 against 5.6e-6 V^2/Hz on a symmetric van der '
                'Pol at 2 f0). Use PAC.pnoise at that frequency.'
                % (int(harmonic), abs(X), _scale))
        if frequency_aware:
            ## `c(f)` per offset -- one bordered adjoint solve each, cached on
            ## `|offset|` so a symmetric sweep pays once per magnitude.  `c(0)`
            ## is `c` exactly, so the near-carrier lineshape is unchanged.
            off = np.asarray(offsets, dtype=float)
            _cache = {}
            flat = []
            for o in np.atleast_1d(off).ravel():
                key = abs(float(o))
                if key not in _cache:
                    _cache[key] = (c if key == 0.0 else
                                   self.frequency_aware_diffusion(pss, key))
                flat.append(float(self.lorentzian(np.array([o]), _cache[key],
                                                  f0, harmonic)[0]))
            Sv = abs(X) ** 2 * np.asarray(flat).reshape(np.shape(off))
        else:
            Sv = abs(X) ** 2 * self.lorentzian(offsets, c, f0, harmonic)
        with np.errstate(divide='ignore'):
            L = 10.0 * np.log10(np.maximum(Sv / max(abs(X) ** 2, 1e-300),
                                           1e-300))
        return Sv, L

    def frequency_aware_diffusion(self, pss, offset):
        """`c(f)` — the phase diffusion constant seen at modulation offset `f`.

        `c(f) = (1/T) integral v_f^H (CY/2) v_f dt` with `v_f` the
        frequency-aware PPV (`PSS.frequency_aware_ppv`, Lai 2008 eq. 23) on
        the equation rows, integrated by the SAME quadrature as
        `diffusion_constant` over the orbit's full grid -- so `c(0)` is `c`
        exactly, not approximately.

        ⚠⚠ WHY IT EXISTS.  The Lorentzian from `c` uses the DC PPV at every
        offset: a noise current is assumed to move the phase instantly.
        Wherever part of that response goes THROUGH a slow mode -- the
        amplitude mode on an orbit with AM-to-PM coupling, or a slow node in
        the source's path -- it is filtered above that mode's corner, and
        the DC-PPV Lorentzian over-states.  With `c(f)` the Lorentzian
        matches pnoise's PM content to ~2 % over 0.3-10 f_amp on van der Pol
        with an asymmetric orbit (the DC PPV: up to 3.3x high) and to 0.4 %
        behind a slow node (the DC PPV: 0.73x and 0.027x); on a symmetric
        orbit `c(f)/c` stays within 1e-3.

        ⚠ Stationary WHITE sources only, like `diffusion_constant`.  Cost:
        one bordered adjoint GMRES per offset (0.25-0.7 s on these fixtures);
        the solve can fail to converge (Lai's own warning about eq. 23), and
        then this raises rather than returning the DC value silently.

        History: `doc/shooting_history.md`, `PAC.frequency_aware_diffusion`.
        """
        self._check_circuit(pss)
        self._refuse_coloured(pss, 'frequency_aware_diffusion')
        self._refuse_driven(pss, 'frequency_aware_diffusion')
        off = abs(float(offset))
        if off == 0.0:
            return self.diffusion_constant(pss)
        _v, fi = pss.frequency_aware_ppv(off)
        base = fi['ppv']
        m = pss.cir.n - 1
        S = np.asarray(fi['samples_eq'])[:, :m]
        h = np.diff(np.asarray(base['times'], dtype=float))
        n = min(len(h), S.shape[0])
        T = float(base['period'])
        cy = 0.5 * np.real(self._cy_reduced(pss, 2.0 * np.pi / float(pss.period)))
        quad = np.real(np.einsum('ij,jk,ik->i', np.conj(S[:n]), cy, S[:n]))
        return float((quad * h[:n]).sum() / T)

    def _warn_above_amplitude_pole(self, offsets, f0):
        """⚠ THE PHASE-ONLY SPECTRUM IS A LOWER BOUND ABOVE `f_amp`.

        `oscillator_spectrum` returns the PHASE contribution only.  A real
        oscillator also carries AMPLITUDE noise, which is suppressed near the
        carrier because the limit cycle restores the amplitude -- but only at
        the amplitude-relaxation rate.  Above the pole where that restoring
        action runs out, amplitude noise stops decaying within a period and
        adds to the total, so this method UNDER-reports (a relayed
        measurement against a commercial simulator's total noise, ~3 dB low
        well above f_amp).

        ⚠⚠ AND THE VALID REGION SHRINKS AS `1/Q`.  With
        `f_amp = -ln(lam2)/(2 pi T)` and `Q = -1/ln(lam2)`,

            f_amp = f0 / (2 pi Q)

        so the better the oscillator, the narrower the band in which its
        phase-only spectrum is the whole answer; at `lam2 = 0.999` it has
        collapsed below ~253 Hz.

        ⚠ THIS IS THE OPPOSITE SIGN FROM THE ERROR `PSS.ppv` ALREADY WARNS
        ABOUT.  That one says the instantaneous phase equation misses slow
        nodes which FILTER device noise, so phase noise is OVER-estimated.
        This one is a second, independent mechanism in which the phase-only
        answer is UNDER-estimated.  Both are live and they are not the same
        effect.

        History: `doc/shooting_history.md`, `PAC._warn_above_amplitude_pole`.
        """
        lam2, certified = getattr(self, '_last_second_multiplier',
                                  (None, None))
        if lam2 is None:
            return
        lam2 = float(lam2)
        ## `lam2 <= 0` is a real or overdamped mode with no relaxation pole to
        ## speak of, and `lam2 >= 1` is not a decaying mode at all -- in both
        ## cases there is no `f_amp` and inventing one would be worse than
        ## silence.
        if not (0.0 < lam2 < 1.0):
            return
        ## `f_amp = -ln(lam2)/(2 pi T)` and `T = 1/f0`.
        f_amp = -np.log(lam2) * float(f0) / (2.0 * np.pi)
        off = np.atleast_1d(np.asarray(offsets, dtype=float))
        worst = float(np.max(np.abs(off))) if off.size else 0.0
        if worst < f_amp:
            return
        warnings.warn(
            'PAC.oscillator_spectrum: this is a PHASE-ONLY spectrum and %g Hz '
            'is above the amplitude-relaxation pole f_amp = %.4g Hz '
            '(lambda_2 = %.6f, f_amp = f0/(2*pi*Q)). Above f_amp the '
            'amplitude noise no longer decays within a period and adds to the '
            'total, so on a half-wave-symmetric orbit the value returned here '
            'is a LOWER BOUND: measured excess of a commercial simulator over '
            'the phase-only prediction is -0.54 dB at 1 kHz and -2.90 dB at '
            '10 kHz for lambda_2 = 0.99. On an ASYMMETRIC orbit it can instead '
            'OVER-state the total (pnoise/phase-only = 0.61 at half-wave '
            'asymmetry 0.10 with frequency_aware=False, confirmed by Monte '
            'Carlo; the default frequency_aware=True corrects the phase part '
            'to ~2 %%) -- use PAC.pnoise for '
            'the total above f_amp. '
            '%sThe valid band scales as 1/Q, so it NARROWS as the oscillator '
            'improves.'
            % (worst, f_amp, lam2,
               ('' if certified is not False else
                'lambda_2 itself is NOT certified here (see '
                "info['second_multiplier_certified']), so f_amp is uncertain "
                'too. ')),
            RuntimeWarning, stacklevel=3)

    @staticmethod
    def am_pm_indices(a, b):
        """Split a sideband pair into AM and PM modulation indices.

        `a` and `b` are the upper and lower sideband amplitudes, each
        already divided by the carrier phasor.  Returns `(m_am, m_pm)`.

        THE WHOLE THING IS ONE CONJUGATE.  Write the complex envelope's
        deviation as `a e^{j w_m t} + b e^{-j w_m t}`.  The two sidebands
        COUNTER-ROTATE about the carrier phasor, so the sum traces an
        ellipse; the component ALONG the carrier is amplitude modulation
        and the component PERPENDICULAR to it is phase modulation.

          - pure AM keeps the envelope on the carrier's axis, which forces
            `d = conj(d)` for all `t`, i.e. `a = conj(b)`;
          - pure PM keeps it perpendicular, `d = -conj(d)`, i.e.
            `a = -conj(b)`.

        so `m_am = a + conj(b)` and `m_pm = a - conj(b)` -- each vanishing
        exactly when the other case holds.  No new solve: this is a change
        of basis on transfer functions `adjoint_sideband_row` already
        returns.

        ⚠ `conj(b)`, NOT `b`.  Using `a +- b` looks equally plausible and
        is wrong for any modulation whose sidebands are not real relative
        to the carrier -- it would report a rotating ellipse as pure AM.
        The conjugate is what makes the lower sideband counter-rotate.
        """
        a = np.asarray(a, dtype=complex)
        b = np.asarray(b, dtype=complex)
        return a + np.conj(b), a - np.conj(b)

    def _output_waveform_row(self, pss, output):
        """The steady-state waveform of `output`, as an index OR a direction.

        An integer names a node, so callers that name one keep working; an
        array is a DIRECTION, contracted against the full waveform with the
        reference row reinserted -- as `pnoise`, `adjoint_transfer_row` and
        `adjoint_sideband_row` accept `d`.  That is what lets `am_pm` and
        `carrier_phasor` express a DIFFERENTIAL output: for an oscillator
        the output of interest is very often differential, and for the
        coordinate-invariance an AM/PM split has to have (Kaertner 1990
        section 3.2) a reference-independent observable is the whole point.

        History: `doc/shooting_history.md`, `PAC._output_waveform_row`.
        """
        if getattr(pss, 'waveform', None) is None:
            raise RuntimeError(
                'PAC: the PSS has no stored waveform -- call solve() first.')
        _times, X = pss.waveform
        Xf = np.asarray(X, dtype=float)
        irn = pss.irefnode
        d = np.asarray(output)
        if d.ndim == 0:
            k = int(d)
            return Xf[k if k < irn else k + 1]
        row = np.zeros(Xf.shape[1], dtype=float)
        for i, wgt in enumerate(np.asarray(d, dtype=float)):
            if wgt != 0.0:
                row = row + wgt * Xf[i if i < irn else i + 1]
        return row

    def carrier_phasor(self, pss, output, carrier=1):
        """The `carrier`-th Fourier coefficient of the steady-state output.

        Computed here rather than taken from `fpss`, whose spectrum is RMS
        and energy-folded -- correct for reporting a magnitude and useless
        for a phasor, since folding discards the phase the AM/PM split is
        made of.
        """
        times, _X = pss.waveform
        row = self._output_waveform_row(pss, output)
        t = np.asarray(times, dtype=float)[:-1]
        v = row[:len(t)]
        w0 = 2.0 * np.pi / float(pss.period)
        ## ⚠ a Fourier INTEGRAL: `1/N` is its quadrature only on a uniform
        ## grid (measured 7.5 % off and not converging on a 3:1 one)
        _wq = pss._period_quadrature(pss.factored_period())
        if _wq is None or len(_wq) != len(t):
            return complex(np.sum(v * np.exp(-1j * carrier * w0 * t)) / len(t))
        return complex(np.sum(v * np.exp(-1j * carrier * w0 * t) * _wq))

    def am_pm(self, pss, freq, output, carrier=1):
        """AM and PM modulation indices at `carrier`, per noise/signal source.

        Returns `(m_am, m_pm)`, each a row of length `m`: the modulation a
        unit source at reduced coordinate `i`, driven at `freq`, imposes on
        the `carrier`-th harmonic of the output.

        ⚠ TWO SOLVES AT ±freq, NOT ONE.  The upper sideband of harmonic `i`
        sits at `i f0 + freq` and the lower at `i f0 - freq`; with the
        convention that an input at `f` produces output at `f + l f0`,
        those are `H_i(freq)` and `H_i(-freq)`.  They are NOT conjugates of
        each other -- that would hold for an LTI circuit, and the whole
        point of an LPTV analysis is that it does not.  Taking one and
        conjugating it would silently force `m_pm = 0` or `m_am = 0`
        depending on which.

        ⚠ AND AN OSCILLATOR IS ALMOST PURE PM NEAR ITS CARRIER, which is
        the physical check: the phase response to a perturbation goes as
        `1/w_m` while the amplitude response stays bounded, so
        `|m_pm|/|m_am|` grows without bound as `freq -> 0`.  A
        decomposition that got the conjugate wrong gives a bounded ratio
        instead.

        ⚠ THE ABSOLUTE MAGNITUDE ON AN OSCILLATOR IS SMALL FOR A REASON.
        These are the p = 0 band of `am_pm_noise`: a source at BASEBAND
        `freq` reaching the carrier sideband.  A baseband current moves the
        PHASE through the PPV's DC coefficient (Hajimiri-Lee's c_0), and a
        half-wave-symmetric orbit -- odd nonlinearity, `u(t + T/2) = -u(t)`
        -- has none, so on such a fixture the rows measure a symmetry zero
        (proportional to 1/freq and to mu), the same zero the coloured
        up-conversion gate records for Gamma; breaking the symmetry lifts
        them linearly in the asymmetry.  The DIRECT rows (source at
        f0 + freq, sideband 0) agree with `pnoise` at every offset.  So do
        not read a small `am_pm` on a symmetric oscillator as a defect: it
        is the 1/f^3 up-conversion coefficient, and it is zero there.

        History: `doc/shooting_history.md`, `PAC.am_pm`.
        """
        C = self.carrier_phasor(pss, output, carrier)
        ## ⚠ RELATIVE TO THE SIGNAL, NOT AGAINST ZERO.  A harmonic the
        ## circuit does not produce still has a phasor of ~1e-16 rather
        ## than exactly 0, and dividing by it turns "there is no carrier
        ## here" into an enormous, confident modulation index.
        row = self._output_waveform_row(pss, output)
        scale = float(np.max(np.abs(row)))
        if abs(C) <= 1e-9 * max(scale, 1e-300):
            raise ValueError(
                'PAC.am_pm: the output carries no component at harmonic %d '
                '(|C| = %.3e against a signal scale of %.3e), so there is '
                'no carrier to modulate and AM/PM are not defined. Dividing '
                'by it would report a huge modulation of nothing. Pick a '
                'harmonic the circuit actually produces.'
                % (carrier, abs(C), scale))
        upper = self.adjoint_sideband_row(pss, freq, output, carrier)[0]
        lower = self.adjoint_sideband_row(pss, -freq, output, carrier)[0]
        return self.am_pm_indices(upper / C, lower / C)

    def am_pm_noise(self, pss, freq, output, carrier=1, maxsidebands=None,
                    modulated=False):
        """Output NOISE split into its AM and PM parts at `freq` from `carrier`.

        Returns `(S_am, S_pm, bands_used)`.  The two add to the noise in the
        pair of sidebands they decompose -- see the identity below -- and are in
        the same units as :meth:`pnoise`.

        ⚠ THIS NEEDS THE SIDEBAND *CORRELATION*, WHICH IS WHY IT IS NOT
        `|m_am|^2` FROM :meth:`am_pm`.  That method is the TRANSFER pair for a
        deterministic input; noise asks a different question, because whether
        the upper and lower sidebands are CORRELATED is exactly what decides the
        split.  Uncorrelated sidebands carry equal AM and PM -- the classical
        result for narrowband noise through an LTI system -- and it is the
        periodic operating point that correlates them.

        THE BAND BOOKKEEPING, which is the whole of the derivation and the one
        place a sign error would produce a plausible wrong answer.
        `adjoint_sideband_row(pss, g, output, l)` is the coefficient at output
        `g + l f0` for a unit source at `g`.  The two output sidebands sit at
        `carrier*f0 ± freq`, so a REAL noise band whose positive-frequency
        component is at `g = freq + p f0` reaches

            the UPPER output at `+g` through sideband `l = carrier - p`,
            the LOWER output at `-g` through sideband `l = carrier + p`,

        the second because a real process has `N(-g) = conj(N(g))` -- and that
        shared realisation IS the correlation.  Both contributions come from ONE
        band, so they are combined coherently; different `p` are different
        bands and are summed in power.  :meth:`am_pm` is exactly the `p = 0`
        term of this sum.

        The split per band is the same conjugate one :meth:`am_pm_indices`
        makes, `a + conj(b)` and `a - conj(b)` -- ⚠ the CONJUGATE, not `a ± b`:
        the sidebands counter-rotate about the carrier, and dropping it reports
        a rotating ellipse as pure AM.

        ⚠ THE GATE IS AN IDENTITY, NOT A TOLERANCE.  `pnoise` at the upper
        sideband folds precisely the bands `g = freq + p f0`, and at the lower
        precisely their negatives, so with the factor of one half below

            S_am + S_pm  ==  pnoise(carrier*f0 + freq) + pnoise(carrier*f0 - freq)

        exactly, because `|a+c|^2 + |a-c|^2 = 2|a|^2 + 2|c|^2` leaves no cross
        term.  A pairing error breaks it, which is what the test asserts.

        ⚠ ON A FREE-RUNNING OSCILLATOR this split sits on the SAME absolute
        scale as `pnoise` (the identity to 1e-12) and as the externally
        certified `oscillator_spectrum` (`S_pm = 4 S_v` at every offset: the
        PM content of the pair IS the Lorentzian, 2 S_v per sideband), with
        `S_am` rising from ~0 below the AM corner `f0/(2 pi Q_lambda)` to
        `S_pm` above it -- the ratio is an exact Lorentzian
        `u^2/(u_c^2 + u^2)` in `u = offset/f0`, `u_c = 1/(2 pi Q_lambda)`,
        `Q_lambda = -1/ln|lambda_2|` -- so the pair total is 4 S_v there and
        8 S_v far out.  "~1e-12 rows" from `am_pm` on a half-wave-symmetric
        fixture are a symmetry zero, see `am_pm`.  Oscillator magnitudes
        from this are trustworthy.

        ⚠ WHAT A MEASUREMENT MUST BE TO BE COMPARED WITH `S_pm`: `S_pm` is
        PM BY QUADRATURE OF THE FUNDAMENTAL'S SIDEBANDS.  A "phase" read by
        a one-period demodulation of the fundamental leaks the other
        harmonics' sidebands through its boxcar (sinc(pi(1 - r)) ~ 0.1 in
        amplitude for the second harmonic's), and a phase read from zero
        crossings converts EVERY harmonic's sidebands; both drift (6 %)
        against this quantity on a harmonic-rich asymmetric orbit.  So
        compare `S_pm` with the fundamental's sideband PM (a spectrum
        analyser's sidebands around f0, or the forward-tone gate
        `test_pnoise_oscillator_pm_matches_a_forward_tone_
        transient_with_no_adjoint`), never with a demodulated or
        crossing-time phase, and BAND WITH BAND (see `band_spread`).

        History: `doc/shooting_history.md`, `PAC.am_pm_noise`.
        """
        self._check_circuit(pss)
        pss = pss._adjoint_host()
        fp = pss.factored_period()
        N = len(fp.steps)
        f0 = 1.0 / float(fp.T)
        lmax = N // 2 if maxsidebands is None else min(int(maxsidebands),
                                                       N // 2)
        cyfn = (self._cy_cycle_averaged if modulated else self._cy_reduced)
        k = int(carrier)
        ## ⚠⚠ THE SPLIT IS TAKEN IN THE CARRIER'S FRAME, NOT THE TIME
        ## ORIGIN'S.  AM is the envelope component ALONG the carrier phasor,
        ## so `a + conj(b)` is right only for a cosine-phased carrier;
        ## unrotated, the split depends on where t = 0 sits, and the
        ## identity cannot see it (`|a_r|`, `|b_r|` are `|a|`, `|b|`).
        ## `am_pm` divides by the COMPLEX carrier phasor instead.  With no
        ## carrier at this harmonic the phase is undefined and the split is
        ## left unrotated, as `am_pm` refuses the same case.
        _C = self.carrier_phasor(pss, output, k)
        _scale = float(np.max(np.abs(self._output_waveform_row(pss, output))))
        _rot = (np.exp(-1j * np.angle(_C))
                if abs(_C) > 1e-9 * max(_scale, 1e-300) else 1.0)
        S_am = 0.0
        S_pm = 0.0
        bands = []
        for p in range(-lmax, lmax + 1):
            g = float(freq) + p * f0
            a = self.adjoint_sideband_row(pss, g, output, k - p)[0]
            b = self.adjoint_sideband_row(pss, -g, output, k + p)[0]
            cy = cyfn(pss, 2.0 * np.pi * g)
            a_r = a * _rot
            b_r = np.conj(b) * np.conj(_rot)
            m_am = a_r + b_r
            m_pm = a_r - b_r
            S_am += 0.5 * float(np.real(m_am @ cy @ np.conj(m_am)))
            S_pm += 0.5 * float(np.real(m_pm @ cy @ np.conj(m_pm)))
            bands.append(p)
        return S_am, S_pm, bands

    ## `|1 - alpha|` below which the deflated answer is left as recovered:
    ## nearer a harmonic than this the plain operator is singular at the
    ## arithmetic (an unstaged radau map's multiplier sits at 1e-11 from 1)
    DEFLATION_REFINE_MIN = 1e-8

    def _deflated_solve(self, pss, alpha, b, transposed=False, tol=None):
        """`(I - alpha M) y = b` on an OSCILLATOR, with the pole taken out.

        ⚠ THE SINGULARITY IS THE ANSWER'S OWN POLE, NOT A NUMERICAL DEFECT,
        and that reframing is what makes the fix obvious.  At `alpha = 1`
        the operator is `I - M`, singular by the unit multiplier, and the
        solution really does diverge -- an oscillator's phase response to a
        perturbation goes as `1/df`.  What is wrong is COMPUTING a `1/eps`
        quantity through a system whose conditioning is also `1/eps`: the
        answer is genuinely large and the digits are genuinely gone.

        So the pole is carried ANALYTICALLY.  With `u`, `v` the right and
        left null vectors of `I - M` -- the orbit tangent and the PPV, both
        of which `ppv()` already returns -- border the system:

            [ I - alpha M   u ] [ w ]   [ b ]
            [     v^T       0 ] [ s ] = [ 0 ]

        Because `v^T (I - alpha M) = (1 - alpha) v^T` and `v^T w = 0`, the
        border variable comes out BOUNDED, `s = (v^T b)/(v^T u)`, with no
        `1/eps` in it.  And since `(I - alpha M) u = (1 - alpha) u`, the
        solution is recovered as

            y = w + s u / (1 - alpha)

        where `1 - alpha = 1 - exp(-j w T)` is the factor that vanishes at
        every harmonic, evaluated in closed form rather than inverted
        numerically.

        The bordered operator's conditioning is FLAT (from 0.3 down to 1e-9
        of `f0` on van der Pol) while the plain one tracks the offset; where
        the plain solve is still trustworthy the two agree, and their
        disagreement grows as `1/df` -- the PLAIN solve losing digits.

        ⚠ IT STILL DIVERGES AT AN EXACT HARMONIC, and it should: `1/(1 -
        alpha)` is then a division by zero, and the physical response is
        unbounded.  What changes is that every offset NEAR a harmonic is
        now well conditioned, which is where phase noise is measured.

        `transposed` solves `(I - alpha M^T) x = b`, whose null space is
        spanned by `v` and whose left null space is spanned by `u`, so the
        borders swap.

        History: `doc/shooting_history.md`, `PAC._deflated_solve`.
        """
        import scipy.sparse.linalg as spla
        fp = pss._state_map()
        n = fp.width
        _v, info = pss.ppv()
        v = np.asarray(_v, dtype=float)
        u = np.asarray(info['tangent_pair'], dtype=float)
        vu = float(v @ u)
        if abs(vu) < 1e-300:
            raise ValueError(
                'PAC: the PPV is orthogonal to the orbit tangent, so the '
                'bordering is singular and the pole cannot be removed.')
        col, row = (v, u) if transposed else (u, v)
        ## ⚠ THE BORDER VECTORS ARE NORMALISED: the recovery
        ## `y = w + s col / (1 - alpha)` is scale-free in exact arithmetic,
        ## but GMRES sees the bordered MATRIX, and with the tangent in V/s
        ## against the PPV in s/V its condition number can be orders above
        ## that of `I - alpha M` alone.  Unit vectors put the bordering at
        ## the operator's own conditioning; `s` absorbs the scale.
        col = col / max(float(np.linalg.norm(col)), 1e-300)
        row = row / max(float(np.linalg.norm(row)), 1e-300)
        mv = (fp.matvec_transposed if transposed else fp.matvec)
        ## ⚠ ON A STAGED SOLVE THE POLE IS THE TOTAL MAP'S: `u`, `v` are the
        ## null vectors of `I - M_tot`, `M_tot = M + P_theta dtheta/dx_0`, and
        ## the fixed-grid `M` has no unit multiplier at all.  The operator
        ## here is the total map's.
        _ev = EventColumns.of(pss, n)
        if _ev is not None:
            mv = _ev.total_matvec(mv, transposed=transposed)
        b = np.asarray(b, dtype=complex).ravel()

        def _mv(z):
            z = np.asarray(z)
            w_, s_ = z[:n], z[n]
            top = w_ - alpha * np.asarray(mv(w_)) + s_ * col
            return np.concatenate((top, [complex(row @ w_)]))

        A = spla.LinearOperator((n + 1, n + 1), matvec=_mv, dtype=complex)
        rhs = np.concatenate((b, [0.0 + 0.0j]))
        rt = max(self.KRYLOV_FACTOR * pss.par.reltol if tol is None else tol,
                 1e-14)
        z = self._gmres_checked(A, rhs, rt, 'the deflated solve')
        w, s = z[:n], z[n]
        denom = 1.0 - alpha
        if denom == 0:
            raise ValueError(
                'PAC: the deflated solve was asked for an EXACT harmonic, '
                'where 1/(1 - alpha) is a division by zero and the physical '
                'response is unbounded. The pole is removed from the '
                'CONDITIONING, not from the answer.')
        y = w + s * col / denom
        ## ⚠ REFINED ON THE PLAIN OPERATOR WHERE THAT IS WELL CONDITIONED:
        ## the recovery assumes `u`, `v` are EXACT null vectors of `I - M`.
        ## On a staged solve the discrete total map's unit multiplier is
        ## displaced by O(h), and the recovered `y` then misses the true
        ## operator by that much over `|1 - alpha|`.  The plain operator is
        ## well conditioned there (its pole sits where the DISCRETE
        ## multiplier is, not at alpha = 1), so the deflated answer seeds a
        ## plain correction on its own residual, kept only if it lowers the
        ## residual; below `DEFLATION_REFINE_MIN` the deflated answer
        ## stands.  ⚠ This makes the forward and adjoint solves the discrete
        ## operator's own, dual-consistent -- and near a harmonic that
        ## answer carries the multiplier's displacement (the staged map's
        ## multiplier, roadmap E8).  On an unstaged oscillator the residual
        ## is already at the tolerance and nothing happens.
        ## ⚠⚠ NOT ON A NORDSIECK GLM'S MAP, MEASURED (2026-09-25; Andreas:
        ## "If GLM is accurate use GLM").  Its map opens with the startup,
        ## which breaks the discrete phase symmetry: the unit multiplier sits
        ## ``eta = O(h^p)`` off 1 (van der Pol in LC form, glm3 2.35e-5 /
        ## 1.1e-6 at 60 / 120 points; radau 3e-11).  Refined, the answer is
        ## the discrete operator's and misses the physical one by ``eta /
        ## (2 pi r)``, `r` the offset in units of f0 (glm3 at 60 points:
        ## 3.8e-3 at 1e-3, 0.35 at 1e-5, and below ``r ~ eta / 2 pi`` the
        ## pole is gone).  Unrefined -- the pole carried analytically -- it
        ## is O(h^p) at every offset (3.3e-5 / 3.3e-5 / 8.4e-5 at 1e-3 /
        ## 1e-5 / 1e-7 against radau at 480 points), and forward and adjoint
        ## then agree to O(eta) rather than the arithmetic (1e-6).
        if abs(denom) >= self.DEFLATION_REFINE_MIN and not fp.is_glm:
            def _plain(z_):
                z_ = np.asarray(z_)
                return z_ - alpha * np.asarray(mv(z_))
            r0 = b - _plain(y)
            scale = max(float(np.linalg.norm(b)), 1e-300)
            if float(np.linalg.norm(r0)) / scale > rt:
                dy, _rr, _H2, _k2 = _arnoldi_gmres(_plain, r0, rtol=rt,
                                                   maxiter=min(n, 200))
                if float(np.linalg.norm(b - _plain(y + dy))) < float(np.linalg.norm(r0)):
                    y = y + dy
        return y

    def _op(self, fp, alpha):
        """`v -> (I - alpha M) v`, never forming `M`."""
        return lambda v: np.asarray(v) - alpha * fp.matvec(v)

    def _solve_each(self, fp, alphas, rhs, tol):
        """One GMRES per frequency -- the baseline the sweep is measured against."""
        import scipy.sparse.linalg as spla
        n = fp.width
        ys, count = [], [0]

        for alpha, b in zip(alphas, rhs):
            op = self._op(fp, alpha)

            def _mv(v, _op=op):
                count[0] += 1
                return _op(v)

            A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
            ys.append(self._gmres_checked(A, b, tol, 'the m x m solve'))
        return ys, count[0]

    def _solve_subspace(self, fp, alphas, rhs, tol):
        """ONE Krylov subspace for the whole sweep -- the recycling.

        ⚠ THE SUBSPACE IS FREQUENCY-INDEPENDENT AND THAT IS THE WHOLE POINT.
        `A(alpha) = I - alpha M`, so

            span{r, A r, A^2 r, ...} = span{r, M r, M^2 r, ...}

        for every `alpha` -- Telichevesky et al.'s Theorem 1.  A basis of
        `M`'s Krylov space therefore serves every frequency, and each one
        costs a small dense least-squares over it instead of its own run of
        full-period replays.

        ⚠ WHAT IS NOT FREE is the right-hand side: `w(f)` genuinely differs
        per frequency, so a basis grown from one frequency's residual is not
        the space GMRES would have chosen for another.  This does not
        guess -- it minimises the TRUE residual over the shared span, checks
        it, and extends the basis (one matvec, kept for every later
        frequency) until every frequency is inside tolerance.  So the answer
        is never worse than the per-frequency solve; only the matvec count
        varies.
        """
        n = fp.width
        V = np.zeros((n, 0), dtype=complex)
        MV = np.zeros((n, 0), dtype=complex)
        count = [0]

        def extend(seed):
            """One Arnoldi step of `M` from `seed`, orthogonal to `V`."""
            nonlocal V, MV
            v = np.asarray(seed, dtype=complex).ravel().copy()
            if V.shape[1]:
                v = v - V @ (V.conj().T @ v)
                v = v - V @ (V.conj().T @ v)   ## reorthogonalise once
            nv = np.linalg.norm(v)
            if nv < 1e-14:
                return False
            v = v / nv
            count[0] += 1
            Mv = np.asarray(fp.matvec(v), dtype=complex)
            V = np.hstack((V, v[:, None]))
            MV = np.hstack((MV, Mv[:, None]))
            return True

        extend(rhs[0])
        ys = [None] * len(alphas)
        pending = list(range(len(alphas)))
        for _round in range(min(n, 200)):
            still = []
            for i in pending:
                alpha, b = alphas[i], rhs[i]
                AV = V - alpha * MV
                y, *_ = np.linalg.lstsq(AV, b, rcond=None)
                x = V @ y
                r = b - (x - alpha * (MV @ y))
                ys[i] = x
                nb = np.linalg.norm(b)
                if np.linalg.norm(r) > tol * (nb if nb else 1.0):
                    still.append((i, r))
            if not still:
                return ys, count[0]
            if V.shape[1] >= n:
                break
            ## grow the shared basis on the worst residual -- one matvec,
            ## and every frequency gets to use it
            worst = max(still, key=lambda p: np.linalg.norm(p[1]))
            if not extend(worst[1]):
                break
            pending = [i for i, _r in still]
        return ys, count[0]
