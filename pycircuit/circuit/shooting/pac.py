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

    THE OPERATOR IS THE MONODROMY, and the whole method is one line of
    algebra on the withdrawn implementation's own system.  That system was

        (L + alpha B) v = -u,   alpha = exp(-2j pi f T)

    with `L` the block lower bidiagonal discretisation over the period and
    `B` the periodic wrap.  Telichevesky, Kundert & White (DAC 1996) reach
    the iterative form by "reinterpreting the use of `L^-1` ... as a
    preconditioner":

        (I + alpha L^-1 B) v = -L^-1 u

    `L` is block lower bidiagonal, so applying `L^-1` is forward
    substitution through the timesteps -- which is the recursion PSS already
    runs against stored factors -- and `B` is confined to the first `m` rows
    and last `m` columns, so `L^-1 B` acts only on the LAST block.  Both
    claims are checked against our own matrices in
    `test_the_pac_operator_is_the_monodromy_and_L_is_never_formed`.

    What is left after that is `m x m`:

        (I - alpha M) y_0 = alpha w(f)

    with `M` the monodromy and `w` the forced response over one period from
    a zero initial state.  `y_0` is the small-signal state at `t = 0`; one
    more driven replay gives the rest of the period.

    ⚠ WHY THE 419.5 GiB IS GONE, precisely.  It was never the operator: it
    was the cost of FORMING `L` and `B`, `(N m)^2` complex entries, 279.7 +
    139.8 GiB at `N = 137`, `m = 1000`.  Nothing here forms either.  The
    stored per-step factors PSS already makes are the preconditioner, and
    the only dense object is `m x m` and only if the caller asks for it.

    ⚠ AND THE OLD `L` WAS BACKWARD-EULER-SHAPED, which is the trap a rewrite
    falls into.  It has two terms per row; a two-step method's variational
    system has three.  Rebuilding it for `trap` or `gear` gives an operator
    for a different recursion than the trajectory it came from -- measured,
    spectral radius 0 against the analytic 0.8546
    (`test_the_pac_L_is_backward_euler_only`).  Taking `M` from the
    traversal cannot make that mistake, because every step carries its own
    `(alphas, b)`.
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
        """
        toolkit = self.toolkit
        freqs = np.atleast_1d(np.asarray(freqs, dtype=float))
        fp = pss.factored_period()
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
        ## takes `epar` second, and the withdrawn body wrote
        ## `self.cir.u(0, analysis_name)` -- passing 'ac' as the element
        ## parameter set and taking the TRANSIENT source vector, which is
        ## zero at `t = 0` for every sinusoid.  The whole analysis would
        ## have returned zeros, silently, with no source to speak of.
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
        ## ORDER.  On the plain path `_traverse_factored_plain` takes one
        ## step OUTSIDE the loop to manufacture a history, and folds its
        ## effect into the `opening` triple as a flat-history assumption.
        ## For the HOMOGENEOUS map that is the documented approximation the
        ## whole plain path is built on.  For the DRIVEN one it also means
        ## the source is never applied at that step -- one step of `u` out
        ## of `N`, i.e. a relative O(h).
        ##
        ## ⚠ MEASURED, on the Q=20 resonator against the AC analysis at
        ## 700 Hz, rel error per doubling of the grid:
        ##
        ##     trap, plain            2.00x  (O(h))   4.13e-03 at 250 pts
        ##     trap, x0_unknown=True  4.00x  (O(h^2)) 1.09e-04 at 250 pts
        ##     euler, either          2.00x  (O(h))   1.40e-02, unchanged
        ##
        ## The euler row is the control: `x0_unknown` does not move it at
        ## all (identical to five digits), so the trapezoidal gain is the
        ## manufacturing step and not something else the formulation does.
        ## The trajectory is NOT the problem -- trap's waveform converges at
        ## 4.2x per doubling either way.
        ##
        ## So this is a silent order loss for a caller who did nothing
        ## wrong, which is the one thing worth a warning.  Gear-2 takes the
        ## solved-history path and has no manufacturing step at all.
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

        rhs = []
        for f in freqs:
            w, _ = pss._forced_replay(fp, f, u_ac)
            rhs.append(np.exp(-2j * np.pi * f * T) * np.asarray(w))

        alphas = [np.exp(-2j * np.pi * f * T) for f in freqs]
        tol = max(pss.par.reltol * self.KRYLOV_FACTOR, 1e-14)
        ## ⚠ ON AN OSCILLATOR THE OPERATOR HAS THE ANSWER'S OWN POLE at every
        ## harmonic (see `_check_harmonic`), and a plain solve near one
        ## carries relative error `eta / (2 pi df/f0)`, `eta = |lambda_1 - 1|`
        ## the computed unit multiplier's displacement -- measured to four
        ## digits over five decades (docs session, Gourary reading).  The
        ## deflated route (`_deflated_solve`) borders the pole out and is
        ## exact there; it was wired into `adjoint_sideband_row` only, and
        ## this sweep solved plain outside HARMONIC_GUARD (2026-09-08).
        ## Under the radau default eta ~ 1e-12 puts the unguarded band
        ## inside the guard, so this is correctness hygiene, not a fix a
        ## user would see; the subspace recycling across frequencies is
        ## given up on the autonomous path (one bordered solve per point).
        self.deflated = bool(getattr(pss, 'autonomous', False))
        _dth_f = [None] * len(freqs)
        if self.deflated:
            ## ⚠ ON A STAGED OSCILLATOR (2026-09-22, events phase B) the
            ## bordered system collapses onto the total map: with `dtheta
            ## = dtheta/dx_0 y_0 + dtheta_f`, `dtheta_f = -Gt^-1 W f_node`
            ## the source's own motion of the crossings, `(I - a M_tot)
            ## y_0 = a (w + P_theta dtheta_f)` -- the deflated solve with
            ## the total operator and this source.  Exact against the
            ## piecewise-linear forced response (see the test); the plain
            ## deflated solve on such a solve was 0.3-400x off.
            _evd = EventColumns.of(pss, fp.width)
            if _evd is not None:
                _Pthd = np.asarray(_evd['P_end'], dtype=complex)
                for i, (f, a) in enumerate(zip(freqs, alphas)):
                    _e0, f_steps = pss._forced_replay(fp, f, u_ac, y0=np.zeros(m, dtype=complex),
                                                      collect=True)
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
        ## ⚠ THE BORDERED SIDEBAND RESPONSE (2026-09-22, events phase B).
        ## On a solve whose grid was landed on state events, a periodic
        ## perturbation moves the crossings: `y_end = M y_0 + w + P_theta
        ## dtheta`, and the event rows close it -- `w_k . y(node_k) = 0`
        ## with `y(node) = P_node y_0 + f_node + Pk_node dtheta` (the
        ## homogeneous map to the node, the forced response there, the
        ## event column there).  Solved by block elimination: the m x m
        ## solve for the source and for each event column, then the K x K
        ## Schur complement for `dtheta`.  A per-step saltation instead of
        ## this read the dominant multiplier 28 % short of the exact total
        ## (see `_state_event_stage`); the bordered system IS the
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
                                                  collect=True)
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
            for i, y0 in enumerate(ys):
                dthetas[i] = _dthx @ np.asarray(y0, dtype=complex)[:m] + _dth_f[i]
        ## the crossings' modulation per frequency (fractions of the period
        ## per unit source), None where the solve had no state events
        self.event_shifts = list(dthetas)
        outfreq, outV = [], []
        ## the complex time-domain response per frequency, `(times, y)` with
        ## `y` the response to the source `u_ac e^{j w t}` at the grid's
        ## nodes (reduced state); the sideband coefficients below are its
        ## periodic envelope's Fourier coefficients by the period quadrature,
        ## which on a strongly non-uniform grid is not an interpolating basis
        ## -- a time-domain reading comes from here, not from summing them
        self.time_response = []
        _Pk_fixed = None
        for f, y0, dth in zip(freqs, ys, dthetas):
            _end, ysteps = pss._forced_replay(fp, f, u_ac, y0=y0, collect=True)
            y = np.array([np.asarray(y0)[:m]] + [np.asarray(v)[:m]
                                                 for v in ysteps])
            if dth is not None:
                ## the crossings' motion at every node, AT FIXED TIME --
                ## `Pk_j - xdot_j tau_j^T` (2026-09-22): the response of "node
                ## j" itself includes the node's motion along the orbit,
                ## O(1) of the response on a staged oscillator; with the
                ## motion removed the exact forced response is matched to
                ## 1e-3 at every node (see `_fixed_time_event_columns`)
                if _Pk_fixed is None:
                    _Pk_fixed, _t_, _x_ = self._fixed_time_event_columns(pss)
                y = y + np.tensordot(_Pk_fixed[:len(y)], dth, axes=(2, 0))
            self.time_response.append((np.asarray(fp.times, dtype=float)[:len(y)], y.copy()))
            ## `v(t) = y(t) exp(-j w t)` is T-periodic; its DFT is the
            ## sideband set, exactly as the withdrawn body intended
            tms = np.asarray(fp.times, dtype=float)[:len(y)]
            v = y * np.exp(-2j * np.pi * f * tms)[:, None]
            ## ⚠ TWO REPORTING DEFECTS, FOUND BY AN EXTERNAL REFERENCE CROSS-CHECK
            ## (2026-09-05), neither in the solve.  (a) `fp.times` spans
            ## `[0, T]` INCLUSIVE, so the last sample repeats the first on a
            ## T-periodic `v` (|v[0] - v[-1]| / |v[0]| = 7e-18 measured) and
            ## the DFT's `dt = T/(N-1)` put the sidebands at `f0 (N-1)/N`:
            ## 99 500 Hz for 100 000 at N = 200 -- and cost an ORDER, O(h)
            ## for O(h^2), 68x at 800 points.  `PSS.solve` already drops
            ## the endpoint one function away; this did not.  Guarded on
            ## the window rather than sliced blind, since the plain path's
            ## `[:len(y)]` need not be inclusive.  (b) `|sb + f|` folded a
            ## NEGATIVE sideband frequency to positive and left the
            ## coefficient alone; the physical response there is the
            ## CONJUGATE.  Uncorrected, `l = -1` was 166% off and did not
            ## converge under refinement; conjugated it lands on its
            ## positive twin's error to three digits (4.873e-3 / 4.877e-3).
            ## Both defects are invisible on a circuit whose `v(t)` is
            ## constant over the period -- every earlier PAC gate.
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
        replay, and the whole row falls out. MEASURED against `m` forward
        solves on an RC ladder: agreement 9.6e-16.

        ⚠ WHAT THIS IS NOT, so the next reader does not over-read it. The
        output here is the state at `t = 0`, a single linear functional.
        A SIDEBAND coefficient `H_l` is a functional DISTRIBUTED over the
        period -- `(1/N) sum_n exp(-j l w0 t_n) d^T y_n` -- and its adjoint
        needs the reverse pass to take an injection at every step rather
        than a seed at the end. That extension is the next piece of A3, and
        it is not built.

        ⚠ WAS SOLVED-HISTORY ONLY until B8 gave the one-step companions
        their own reverse recursion; it now runs under every method.
        """
        import scipy.sparse.linalg as spla
        ## ⚠ NO LONGER GEAR-ONLY (B8): the transposed replay exists for the
        ## one-step companions too, and every use below goes through
        ## `fp.matvec_transposed`.
        fp = pss.factored_period()
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
        Dropping the second term would leave an answer that looks entirely
        reasonable: MEASURED on an RC ladder the two terms are comparable
        in size (303 against 498 at `l = 0`, 79 against 606 at `l = 1`), so
        neither is a correction to the other.

        Still ONE transposed solve per sideband whatever the number of
        sources, which is the property pnoise needs.  Agreement with the
        `m` forward driven solves: 9.2e-16 / 3.3e-16 / 1.3e-15 at
        `l = 0 / 1 / -2`.

        ⚠ WAS SOLVED-HISTORY ONLY, like the reverse pass; B8 lifted both.
        """
        import scipy.sparse.linalg as spla
        ## ⚠ NO LONGER GEAR-ONLY (B8) -- see the adjoint row.
        fp = pss.factored_period()

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
            ## ⚠ ON A STAGED SOLVE THE ROW IS BORDERED (events phase B,
            ## 2026-09-22): the transpose of `PAC.solve`'s bordered system,
            ## the output read at FIXED times (`g_theta` over the fixed-time
            ## columns), and the event rows' term -- `-zeta_k W_k` at node
            ## k, the source coupling of `f_node_k` -- as a second reverse
            ## pass.  On a driven solve `z, zeta` come from the block
            ## elimination (`EventColumns.bordered_adjoint`); on an
            ## oscillator the system collapses onto the TOTAL operator,
            ## deflated, with `zeta` read off after it.  Radau/trbdf2 and
            ## gear's pair map alike (an autonomous gear solve stays
            ## unbordered).  Verified by dual consistency against the
            ## bordered forward solve (the suite's own PAC test pattern);
            ## unbordered, pnoise on a staged gear solve was 10-15 % off.
            _autonomous = getattr(pss, 'autonomous', False)
            _ev = EventColumns.of(pss)
            if _ev is not None and not (fp.is_stage
                                        or (fp.is_pair and not _autonomous)):
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
    ## `oscillator_spectrum` for why the two cannot be unified and why the wrong
    def pnoise(self, pss, freq, output, ratio_tol=None, maxsidebands=None,
               modulated=False, cyclostationary=False):
        """TIME-AVERAGED output noise PSD at `freq`, sidebands folded in.

        ⚠ `cyclostationary=True` IS THE CONSTRUCTION FOR A BIAS-DEPENDENT
        `CY` (2026-09-08, on the corrected Okumura reading).  A source whose
        PSD follows the orbit is white noise `xi` MODULATED by
        `B(t) = sqrt(CY(x(t)))`, a T-periodic matrix with Fourier
        coefficients `B_k` -- read off the PSS samples by one DFT, no
        window count `p` at all (Okumura's windows are a piecewise-constant
        approximation of exactly this, and their boxcar coefficients its
        crude version).  `xi`'s band at `g_p = f - p f0` reaches the output
        at `f` through EVERY modulation harmonic `k` and the sideband row
        `a_{p-k}` (source at `g_p + k f0`, output at `f`) -- the SAME rows
        the stationary fold computes -- COHERENTLY over `k` (one white band,
        one realisation) and incoherently over `p`.  Summing the bands
        turns the square root into the PSD's OWN harmonics `P_j` (the DFT
        of `CY(x(t))`, no matrix square root anywhere):

            S(f) = sum_{l,l'} a_l P_{l'-l} a_{l'}^H,

        which is exact on the grid (⚠ the sqrt-modulation form, tried
        first, left a 2.8e-5 residual tied to the modulation's zero
        crossings; this form agrees with the stationary side to 9e-16).
        Constant `CY` gives `P_0 = CY` and nothing else, and the sum
        collapses to the stationary `sum_l a_l CY a_l^H` -- Okumura's
        `p = 1` case, pinned to machine precision.  The cost is the
        stationary fold's (the rows dominate; the double sum is free).
        A coloured source is folded band by band (each white band the
        modulation reaches carries its own `CY`; see `_cyclostationary_fold`).
        Like the stationary fold this is a LOWER bound at a sideband cap.
        ⚠ Coherence is the whole content: `modulated=True` (the cycle-
        averaged `CY`, Hull & Meyer's stationary equivalent) keeps the
        power and drops the correlation between sidebands, and the two
        differ wherever the modulation has harmonics -- measured on a
        driven multiplier, and the identity against the STATIONARY fold of
        the same physics written as a white source through a periodically
        varying gain is the gate (`test_..._cyclostationary_...`).  ON A
        MOS (2026-09-09): an EKV stage switched by a 1 MHz LO whose channel
        noise passes through a second EKV switched by the same LO reads
        0.376 of the cycle average (thermal, every offset) and 0.32-0.44
        with flicker at ten times thermal; the switch's own channel noise
        is largest when its channel shunts it.  ⚠ Noise that reaches the
        output through a time-INVARIANT transfer (a single stage's drain
        into an RC load) gives cyc = cycle average to four digits, since
        only P_0 survives -- the construction shows only where the
        modulated noise crosses a periodically varying transfer.  Cost:
        white = the cycle average's; coloured (any frequency-dependent
        CY, a negligible flicker coefficient included) ~6x.
        ⚠ FLICKER, AND WHAT OKUMURA'S EQ. 23 MEANS HERE (measured
        2026-09-08): a coloured source is folded band by band, and against
        the stationary fold of the same SEPARABLE physics (a stationary
        flicker source through a periodically varying gain) it is exact --
        1.000000 -- as long as the modulation is SIGN-DEFINITE.  When the
        modulation changes sign the two are DIFFERENT physics (0.56 / 1.33,
        grid-independent to six digits): for white noise `m xi` and `|m| xi`
        are one process, for a coloured one whose correlation spans the
        sign change they are not (`R(t,t') = m(t) m(t') R_c(t-t')` keeps
        the sign product), and a PSD cannot carry the sign -- so this fold,
        like the HDL model feeding it, is the `|m|` one.  That is Okumura's
        "cannot be modeled as a cyclostationary process by using this
        method, because it has very long time constants" in concrete form.
        A flicker source with a bias-dependent coefficient gets the `|m|`
        number, correct when its modulation does not change sign.

        Returns `(S, sidebands_used)`.  `S` is the one-sided
        **time-averaged** PSD at the output, in the same units as
        `analysis_ss.Noise`'s `Svnout`.

        ⚠ "TIME-AVERAGED" IS NOT A HEDGE, IT IS THE SPECIFICATION, and
        saying so is the whole of this paragraph's job.  TWO SEPARATE
        MECHANISMS make an output noise cyclostationary, and only one of
        them is about the sources:

          1. bias-dependent sources modulated by the time-varying operating
             point -- this is what `_cy_reduced` refuses, because the
             stationary sum would be the wrong model;
          2. the PERIODIC SOURCE-TO-OUTPUT TRANSFER FUNCTION -- which
             applies even when every source is stationary.  A circuit whose
             only noise is the thermal noise of constant resistors STILL
             has cyclostationary output noise.

        The sideband sum here handles (2) correctly and returns its TIME
        AVERAGE.  That is the right answer for most uses and it is
        incomplete for two ordinary RF topologies, both named by Kundert
        (*Introduction to RF Simulation*, v2 2003 -- relayed from the docs
        session, cited not verified here): a NONLINEAR SUBSEQUENT STAGE
        ("an oscillator drives a limiter ... the same is true when an
        oscillator drives a mixer"), and CASCADED STAGES OFF A SHARED
        REFERENCE, where "the second mixer is synchronous with, and tracks
        the variations in, the cyclostationary noise of the first."  The
        test is whether anything downstream can track the PSD's variation:
        if it cannot, the phase is unknown to it and the time average is
        sufficient.

        ⚠ AND A SCALAR CANNOT CARRY WHAT IS MISSING.  Cyclostationary noise
        is CORRELATED between frequencies separated by `k f0`, where
        stationary noise has no correlation between different frequencies
        at all.  This returns one number per output frequency, so it does
        not represent that correlation -- deliberately, and stated here
        rather than left for a caller to discover by getting a wrong answer
        in one of the two topologies above.

            S(f) = sum_l  h_l CY h_l^H ,   h_l = H_l(f - l f0)

        Noise entering at `f - l f0` leaves at `f` through sideband `l`, and
        white sources in disjoint bands are uncorrelated, so the bands add
        in POWER.  Each `h_l` is one adjoint row -- one transposed solve for
        every source in the circuit, which is the whole reason this is
        affordable.

        ⚠ A PRECONDITION FOR THE FIRST COLOURED SOURCE, recorded here
        because it is unreachable today and will be silent when it is not.
        A 1/f source is singular at DC, and folding puts a copy of that
        singularity at EVERY harmonic.  A commercial RF simulator: "place a cluster of
        frequencies near each harmonic ... but AVOID PUTTING FREQUENCY
        POINTS PRECISELY ON THE HARMONICS ... you run the risk of
        generating absurd noise totals because a very narrow noise peak
        artificially has its apparent width greatly magnified by a large
        frequency, and has its amplitude exaggerated by placing a point
        precisely at the singularity."  Plausible nonsense, no error
        raised.  Every source in the discrete library is white, so `freq`
        landing on `k f0` is harmless now; it stops being harmless the day
        one is not.

        ⚠ AND AN OSCILLATOR IS NOT THIS FUNCTION'S PROBLEM AT ALL.  A
        driven circuit's output noise is cyclostationary; an AUTONOMOUS
        one's is STATIONARY, and structurally so -- "cyclostationarity in
        the oscillator's output would, by definition, imply a time
        reference ... noisy autonomous systems cannot provide a perfect
        time reference" (Demir 2002).  That is the physical counterpart of
        `I - M kron M` being exactly singular for an oscillator
        (`test_no_periodic_covariance_exists_for_an_oscillator`): there is
        no cyclostationary object to compute, not a hard one.  Oscillator
        phase noise is a different output shape entirely -- a closed form
        in a few scalars with no frequency sweep -- and is not built.

        ⚠ WHY CYCLOSTATIONARY IS NOT BUILT -- AND THE REASON RECORDED HERE
        FIRST WAS WRONG.  This said the cross terms need "the `R_{m,n}`
        construction from section III-B", unread, as though the window
        Fourier coefficients were an exotic object.  They are not.  The
        model is `c(t) = sum_m n_m(t) w_m(t)` with `w_m` a T-periodic
        RECTANGULAR window over interval `m`, non-overlapping -- so
        `W_{m,k}` is the Fourier series of a BOXCAR, closed form, a `sinc`
        times a phase.  The `n_m` are taken UNCORRELATED across intervals,
        justified because `H(jw,t)` is time-invariant within each one, so
        the sum is INCOHERENT over `m` and coherent only over `k` within a
        single interval.  Nothing there is missing.
        ⚠ THE ACTUAL BARRIER IS COST, which is a different decision -- and
        the cost as first recorded here was OVER-STATED (verified at the
        source by the docs session, 2026-09-08, Okumura et al. 1993).  The
        source count is p x (noisy devices) where p is the number of
        intervals over which "H(jw,t) is time-invariant within each
        interval" -- set by how fast the transfer varies (their Fig. 2 has
        FIVE windows; a switching circuit moves fast only at transitions),
        NOT by the integration grid: the earlier "500-point grid x 50
        devices = 25 000 sources" tied p to the timestep and was high by
        (timepoints)/p, an order or two.  The reported noise analysis ran at
        14.1x the PSS per frequency point (1086 s vs 77 s), "because all
        aliasing components need to be computed" -- and the NEXT sentence,
        elided before: "it is expected that this problem can be greatly
        alleviated using a vectorization technique, because most of the
        computational power is used to solve linear problems" -- which is
        the batched JAX path this tree already carries.  Whether that is
        affordable is unmeasured here.
        ⚠ AND THE METHOD CANNOT MODEL FLICKER (p. 585, verbatim): "Flicker
        noise generated under a periodic large signal excitation cannot be
        modeled as a cyclostationary process by using this method, because
        it has very long time constants and thus equation (23) does not
        hold" -- eq. 23 being the uncorrelated-across-intervals assumption.
        Their fallback is that flicker "may exist as independent noise
        sources which are practically modeled as stationary random
        processes".  So the construction covers cyclostationary thermal
        and shot noise, and NOT one of the three mechanisms `_cy_reduced`
        names it as the precondition for.
        ⚠ AND ITS AUTHORS LEFT THE PHYSICS OPEN: "it is further necessary
        to discuss the correspondence between the actual physical phenomena
        of noises and this modeling".  The windowed-stationary
        decomposition is a numerical construct, and its fidelity to a real
        device is not settled by its numerical validation.
        (Verified at the source 2026-09-08.  The boxcar is the paper's own
        closed form, R_{m,n} = (h_m/T) Sa(n w_s h_m/2) exp(-j n w_s (tau_{m-1}
        + h_m/2)), p. 585 -- an earlier line here claimed the observation as
        ours.)

        ⚠ STATIONARY SOURCES ONLY, AND IT CHECKS -- mechanism (1) above.
        Okumura's cyclostationary model windows each source to a single
        timestep, and the windows'
        Fourier coefficients then CORRELATE the sidebands -- they stop
        adding in power, and the cross terms need the `R_{m,n}` construction
        from his §III-B.  Every noise source in this element library is
        bias-INdependent (a resistor's `4kT/R` does not read `x` at all), so
        the stationary formula is exact for them; a compact device with a
        bias-dependent `CY` is not covered, and this raises rather than
        returning a number that is quietly the wrong model.

        ⚠ `maxsidebands` IS AN ACCURACY KNOB HERE AND A REPORTING KNOB IN
        `PAC.solve`, WHICH IS THE OPPOSITE OF HOW IT READS.  A commercial RF simulator's
        own documentation states the inversion (relayed, cited not verified
        here): reducing sidebands "affects only the amount of information
        generated, not its quality.  HOWEVER, NOISE SOURCES GENERATE
        SIGNALS AT ALL FREQUENCIES, and therefore with PNoise, reducing the
        number of sidebands acts to REDUCE THE NUMBER OF NOISE
        CONTRIBUTIONS in the output and so REDUCES THE ACCURACY of the
        result."  A driven signal lives at the frequencies it is driven at,
        so dropping sidebands drops answers you did not ask for; noise
        lives at all of them, so dropping sidebands drops power that
        belonged in the total.  Capping it here always makes `S` a LOWER
        bound, never a cheaper estimate of the same number.

        ⚠ AN OBSERVABLE SYMPTOM WORTH KNOWING BEFORE IT IS SEEN.  For an
        oscillator `Phi(T) - I` is singular and its null vector IS THE PPV,
        so a near-carrier noise computation is ill-conditioned by
        construction.  Gourary et al. name what that looks like: "the
        standard time domain noise analysis yields FLAT PSD CURVES OR
        CURVES WITH UNEXPECTED SLOPE NEAR THE OSCILLATION FREQUENCY."  If
        oscillator noise ever comes out flat near the carrier, that is the
        singularity -- not the physics, the noise models or the source
        definitions -- which points at the right layer immediately.  The
        published removal (Gourary et al., eq. 27/28: replace the output
        row of J^T by u^T; verified at the source by the docs session,
        2026-09-08) IS built here as `_deflated_solve`, which borders with
        BOTH null vectors and is the better conditioned of the two; it is
        wired into `adjoint_sideband_row` (so into this method) and, since
        2026-09-08, into `PAC.solve` and `adjoint_transfer_row` as well
        (`PAC.deflated` says which route ran).  The plain solve's relative
        error near a harmonic is `eta / (2 pi df/f0)` with `eta =
        |lambda_1 - 1|` the computed unit multiplier's displacement --
        measured 1.1e-12 (Q = 16) and 1.8e-13 (Q = 100) under radau, so
        under the default integrator the unguarded band sat inside
        `HARMONIC_GUARD` and the wiring is hygiene; under gear at
        df/f0 = 1e-10 the plain solve refuses outright (GMRES residual
        1.7e-6) where the deflated one carries the pole to 1 %.

        ⚠ TWO STOPPING RULES, AND THE BOUND IS NOT THE RATIO TEST.  The
        accumulation stops when a sideband pair adds less than `ratio_tol`
        of the running total -- and it can never pass `|l| <= N/2`, the
        grid's own Nyquist, because nothing aliases down from above the
        maximum frequency the grid represents (eq. 32).  An implementation
        with only the ratio test terminates for the wrong reason and, on a
        coarse grid, after summing harmonics its own grid cannot carry.

        GATED against `analysis_ss.Noise` on a linear circuit, where the
        sidebands vanish and this must reduce to the stationary answer --
        Okumura's own `p = 1` case, "exactly the same as that derived for a
        stationary noise".  Measured ratio 1.000000, with every `l != 0`
        contributing ~1e-32 of the total.

        ⚠ AND THE TIME-AVERAGE CHOICE MATCHES THE REFERENCE IMPLEMENTATION,
        which is worth recording because it was documented above as a
        deliberate scope decision and could have been the wrong one.
        a commercial RF simulator's theory notes on PNoise and QPnoise, both: "THE TIME-AVERAGE of
        the noise at the output of the circuit is computed in the form of a
        spectral density versus frequency."  Same quantity, same
        limitation.  (Relayed from the docs session; cited, not verified
        here.)
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
            ## itself is the convolution after the rows are gathered.
            ## The colour model (see `_cy_colour_model`) is fitted ONCE
            ## here and serves both: the 34 orbit sweeps of the stop rule
            ## were 1.1 s of a 4.1 s call after the fold itself was cut.
            ## per-element components, so independent sources ADD in the
            ## coloured fold (see `_cy_components_model`); its call is the
            ## summed model the stop rule reads, as before
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
        ## ZERO.  A `1/f` term is infinite there; and MEASURED on
        ## `PspMosLongChannel`, a flicker term with its coefficient set to
        ## ZERO is `0/0` and returns `nan`:
        ##
        ##     fnt=1, nfa=0        CY(f=0) = nan   <- DISABLED flicker
        ##     fnt=1, nfa=8e22     CY(f=0) = inf   <- the real singularity
        ##
        ## ⚠ THE FIRST IS THE NASTIER ONE: a caller who sets `nfa = 0`
        ## believing flicker is off still gets `nan` out of `pnoise`, with
        ## no exception anywhere.
        ##
        ## ⚠ AND IT IS NOT "HARMONICS ARE BAD".  A driven divider with
        ## white sources returns 1.490351e-17 at exactly `f0`, and at
        ## `2 f0` and `3 f0` -- the fold to DC is harmless when the
        ## sources are defined there.  So this checks the SOURCES at the
        ## frequency that will actually be used, rather than refusing a
        ## harmonic on principle.
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
            ## ⚠⚠ A FINITE PROBE IS NOT A SAFE ONE (peer report, 2026-09-15).
            ## `1/T` rounds (99999.999999999985 Hz for T = 1e-5 s), so at
            ## `f = f0` the folded band sits 1.5e-11 Hz from DC, not ON it: a
            ## 1/f source there is finite and enormous, and the cyclostationary
            ## fold returned 6.3e-2 V^2/Hz against 9.2e-15 at 0.1 % either
            ## side, with no warning.  So a frequency-DEPENDENT CY at the folded
            ## band refuses too; a white one stays allowed (its fold to DC is
            ## harmless, the divider value above).
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

        ## ⚠ AND THE STEEP REGION BESIDE IT IS A SWEEP HAZARD RATHER THAN
        ## A WRONG NUMBER, so it warns instead of raising.  MEASURED with
        ## a real flicker source: the plateau is 1.321483e-16, `f0 + 1` Hz
        ## gives 1.321766e-16 and `f0 + 0.01` Hz gives 1.350047e-16 -- 2%
        ## high, finite, entirely plausible.  The VALUE is right; a grid
        ## that lands there by accident integrates a spike it never
        ## resolved.  A commercial RF simulator: "you run the risk of generating absurd
        ## noise totals because a very narrow noise peak artificially has
        ## its apparent width greatly magnified".
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
        ## A strongly switching circuit does this readily: measured on a
        ## driven diode at 80 points per period, the accumulation reached
        ## l = +-39 without the ratio test ever firing, while folding was
        ## already contributing 62% of the total.
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
        PSD's own harmonics (`a P a^H`), which is exact on the grid, where
        a sqrt-modulation route (tried first) left a 2.8e-5 residual: the
        square root of a PSD that crosses zero has a kink, its harmonic
        tail decays slowly, and the convolution's window -- the sidebands
        the ratio stop kept, 7 here -- truncated it (measured -1.6e-4 /
        -2.8e-5 / -3e-7 at 5 / 7 / 17 sidebands).  The PSD's harmonics
        decay fast, so this form is exact at any window."""
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
            ## made this NaN and silently disabled pnoise's ratio stop
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
        """⚠ A FALLBACK THAT REPRODUCES THE OLD ANSWER LOOKS LIKE AGREEMENT
        (peer, 2026-09-19: found only by a bit-for-bit A/B against an older
        commit).  If an element STATED its signed amplitudes and the fold is
        about to factor that component by sqrt(PSD) anyway, say so."""
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

        ⚠⚠ WHY (2026-09-15, measured).  The coloured fold took ONE symmetric
        square root of the SUMMED `CY(x(t), w)` per band.  Independent
        sources whose modulations differ do not add under a joint root:
        `sqrt(A(t) + B)` cross-couples them at `(t, t')`.  Switch white
        noise `4kT g(t)` plus a constant 1/f source at the same node read
        pnoise(both) = white + flicker + 7.3 % of the total at 0.013 f0
        (+2.2 % at 0.137 f0), and the sampled variance in hold +9.3 % with
        white and flicker inside ONE element.  One root per component
        restores additivity: 1.7e-16 (separate elements) and 2.3e-11 (one
        element, split) against the sum.

        Returns a callable `w -> (N, m, m)` (the summed `CY`, the contract
        `_cy_colour_model` had) carrying `white` (summed), `white_parts`
        `[(key, A)]`, `flicker` `[(key, B, EF)]` (`C = A + B (w1/w)^EF` per
        element, fitted and verified as `_colour_fit`), `perband` (keys
        whose colour did not fit: evaluated per band) and `w1`.  ⚠ Within
        one element all white terms share a root, as do all power-law terms
        -- independence is resolved to element x {white, coloured}.
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
            ## (f - l f0 < 0), and a non-integer power of a negative base is
            ## NaN -- the ratio stop then never fired and every coloured call
            ## ran to the Nyquist bound with a spurious warning (2026-09-15)
            ## ⚠ and numpy division: exactly ON a harmonic the folded band is
            ## w = 0, where a Python float division raised ZeroDivisionError
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
        ## ⚠ THE SIGN (2026-09-19): where the element states its coloured
        ## AMPLITUDES, the folds factor the component with them instead of
        ## with sqrt(PSD) -- see `Element.noise_amplitudes` (hdl.py)
        model.amplitude = self._signed_amplitudes(pss, w1, flicker, states)
        model.w1 = w1
        return model

    @staticmethod
    def _uniform_exponent(B, EF):
        """The one power-law exponent of a component, or None when its
        non-zero entries carry different ones (then `sqrt(B (w1/w)^EF)` is
        not `(w1/w)^(ef/2) sqrt(B)` and must be taken per band)."""
        ## ⚠ ONLY ENTRIES THAT CARRY WEIGHT VOTE (2026-09-19).  The exponent is
        ## fitted from differences of `CY`, so an entry at 1e-12 of the
        ## component's scale -- a MOS flicker source at the sample where Vds
        ## crosses zero -- has its exponent in the rounding of the white part
        ## beside it (measured 1 - 2.2e-09 on an entry of 1.4e-31 against
        ## 7.2e-20).  That one entry failed the whole component into the
        ## per-band route at ONE clock amplitude of a sweep.
        ## ⚠⚠ AND A WEIGHT CUT-OFF WAS THE WRONG REPAIR (same day, peer, 1000
        ## points): my first fix let entries above 1e-9 of the scale vote, and
        ## a finer orbit landed a sample at 7.1e-09 with its exponent off by
        ## 8.1e-09 -- the fallback fired again at two amplitudes and the new
        ## commit reproduced the OLD one to the last bit there.  The exponent's
        ## noise goes as 1/weight, so no cut-off separates them.  What matters
        ## is what a wrong exponent COSTS: giving entry i the exponent `ref`
        ## misstates the component by `r_i |(w1/w)^d_i - 1| ~ r_i d_i |ln(w1/w)|`
        ## of its scale (`r` relative weight, `d` deviation).  Bounded over 50
        ## e-folds of band frequency and held to 1e-9: a genuinely different
        ## exponent (d ~ 1) still fails from a weight of 2e-11 up, while the
        ## two measured offenders cost 2e-19 and 3e-15.  The reference is the
        ## LARGEST entry's exponent, not the first's.
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
        `CY(x(t))` itself -- no square root, exact on the grid (9e-16
        against the stationary fold of the same physics).  COLOURED source
        (`CY` depends on `w`; detected by comparing two bands): the white
        band `p = l + k` shared by rows `l` and `l'` carries its OWN `CY`,
        so `Q_{l,l'} = sum_k B_k^{(l+k)} B_{k+l-l'}^{(l+k) H}` with
        `B^{(p)}` the sqrt-DFT at the band's frequency `|f - p f0|`, summed
        over ALL `N` modulation harmonics `k` (which is what makes the
        square root exact here: the 2.8e-5 of the first version came from
        a window on `k`, not from the root).  ⚠ Measured by the docs
        session on a flicker source: `||P_0||` differs 24x across the bands
        the fold sums, so "the band of l" (the first version's shortcut)
        was a 24x approximation on the case the feature exists for; the
        band-resolved form is pinned against the stationary fold of a
        stationary FLICKER source through the same multiplier.

        COST (2026-09-09): the coloured branch was 6x the white one
        because of the circuit's `CY` (231 bands x 230 samples), not the
        algebra.  With the colour model (`_cy_colour_model`, fitted once
        in `pnoise` and shared with the stop rule) and the pair sum
        vectorised it is 2.2x the white call and below the cycle average
        (1.7 s / 0.8 s / 2.0 s on the switched EKV fixture), exact to
        1e-11 against the per-band evaluation, which remains the fallback
        for a colour the model does not fit."""
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
        ## ⚠ THE COST WAS THE CIRCUIT'S CY, NOT THE ALGEBRA (profiled
        ## 2026-09-09 on a switched EKV stage: 231 bands x 230 samples =
        ## 53 000 CY evaluations, 7.8 s of a 10.4 s fold; the eigen-
        ## decompositions 0.85 s).  Every colour in the library is thermal
        ## plus flicker in 1/f^ef, so THREE evaluations per sample fix each
        ## entry's shape (A + B (w1/w)^ef, ef by a root find on the ratio of
        ## differences), a FOURTH frequency verifies the fit to 1e-8, and
        ## all the bands come from the model with no further circuit
        ## calls; a source whose colour is not of that shape fails the
        ## check and gets the full evaluation as before.
        Nn = Pa.shape[0]
        if model is None:
            model = self._cy_components_model(pss, f, f0)
        ## ⚠ A SPECIFICATION LIMIT, NOT AN IMPLEMENTATION ONE (docs session,
        ## 2026-09-08): a coloured source under a modulation that CHANGES
        ## SIGN is not representable by any fold built from a PSD -- the
        ## correlation R(t,t') = m(t) m(t') R_c(t-t') keeps the sign product
        ## and CY cannot carry it -- so this fold, like the HDL model
        ## feeding it, computes the |m| process (measured 0.56 / 1.33 of
        ## the signed one on a flicker source through a zero-crossing
        ## gain, grid-independent; 1.000000000 for a sign-definite gain).
        ## The sign is invisible here; its NECESSARY condition is a PSD
        ## that touches zero along the orbit with a KINK in its square
        ## root, so that is warned on.  ⚠⚠ THE 2026-09-09 SCOPE NOTE THAT
        ## STOOD HERE WAS WRONG AND IS WITHDRAWN (2026-09-19): it said a
        ## DEVICE's own flicker "has no sign to lose -- sqrt(PSD(x(t))) >= 0
        ## IS the process".  A 1/f current is a slow conductance fluctuation
        ## TIMES the current and follows its sign; on a PSP sampler whose Vds
        ## crosses zero while it conducts that is +0.1 % at one clock
        ## amplitude and 400x at another, against a commercial simulator.
        ## Where the element states its signed amplitudes
        ## (`Element.noise_amplitudes`) the folds use them and nothing below
        ## applies; this limit is for sources that state none.  (White
        ## sources are untouched: uncorrelated across the period, no sign
        ## product survives.)  Okumura's eq. 23
        ## objection to flicker is then the separate, physical question of
        ## whether a trap process is "modulated coloured noise" at all.
        ## The proxy's threshold: a zero crossing SAMPLED on an N-point grid
        ## bottoms out near (pi/N)^2 of the maximum (6e-4 at 200 points on
        ## the gate fixture), while a sign-definite PSD with a ten-fold
        ## swing sits at 1e-2 -- so 1e-2 separates them here; a heuristic,
        ## and it is a warning for that reason.
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
        ## ⚠ THE ORDER OF THE ZERO (peer): a LINEAR sign crossing m ~ a t
        ## gives sqrt(PSD) ~ |a t|, a first-derivative KINK; a sign-definite
        ## quadratic touch m ~ b t^2 gives sqrt(PSD) ~ b t^2, SMOOTH.  The
        ## circular second difference of sqrt(PSD) is 2|a|h at a kink and
        ## 2b h^2 where smooth -- both shrink under refinement, the smooth
        ## one faster -- so a RAW threshold encodes the grid (5e-3 was safe
        ## at 240 points and a false positive below ~100; peer).  Divided
        ## by h/T and by the maximum it is a DERIVATIVE JUMP, grid-
        ## independent at a kink (2|a|T/s_max ~ 4 pi for a sinusoidal
        ## slope, 12.6 here) and falling as h/T where smooth (~2 (2 pi)^2
        ## h/T: 0.33 at 240 points, 1.0 at 80, 2.0 at 40), so 3 separates
        ## them down to ~50 points per period and the separation grows
        ## with refinement.  Clears the squared-gain case (k V_lo^2, exact
        ## to nine digits) that the touch test alone flagged.  ⚠ STILL
        ## NECESSARY, NOT SUFFICIENT, AND THE DETECTOR'S SENSITIVITY RUNS
        ## INVERSE TO THE EFFECT (peer): the indicator is 12.57 for a
        ## sinusoidal crossing, 0.24 for sign|sin|^1.5 and 0 for sign|sin|^2
        ## -- all sign-changing -- while the discrepancy stays O(1):
        ## MEASURED on the flicker identity with the LO shaped to
        ## v |v|^(p-1), B/A = 0.187 / 1.895 (p = 1, warned), 0.204 / 1.779
        ## (p = 1.5, silent), 0.217 / 1.699 (p = 2, silent) at 0.13 / 1.37
        ## f0; steeper crossings (p = 0.5: 0.161 / 2.046, p = 0.8: 0.175 /
        ## 1.970) err MORE and are caught.  So detector and effect are
        ## aligned for p <= 1 and the silent region is exactly p > 1 (the
        ## crossing flatter than linear): the deviation stays O(1) there
        ## while the indicator falls by orders.  A quiet warning is
        ## therefore not evidence of a small discrepancy.
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
        ## ⚠⚠ ONE SQUARE ROOT PER INDEPENDENT COMPONENT, NOT OF THE SUM
        ## (2026-09-15).  A joint `sqrt(CY)` made independent sources with
        ## different modulations NON-ADDITIVE (+7.3 % of the total on a
        ## switch + a 1/f source at one node) -- see `_cy_components_model`.
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
            ## Python (204 000 `_B` calls, 1.9 s of a 2.8 s fold, before).
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
                ## k + l - l', NOT k + l' - l -- the mirror was invisible to
                ## the constant-modulation reduction (only k = 0 there) and
                ## read 0.49 / 0.17 on the smooth-modulation flicker identity.
                ## ⚠ NO CIRCULAR WRAP HERE: a partner beyond N/2 would be
                ## paired with the wrong BAND (each band carries its own
                ## weight), harmless in the white P-form and wrong here --
                ## it read 0.56 / 1.33 on the kinked (zero-crossing)
                ## modulation whose coefficients reach N/2.
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
        PERCENTAGE.  Measured by an external reference-simulator cross-check (2026-09-05)
        on a series switch + shunt capacitor, `pnoise` at 10 kHz against
        a reference simulator, swept over the modulation depth `goff/gon`:

            goff/gon   1        1e-1     1e-2     1e-3     1e-6
            ratio      1.000    4.33     13.2     15.7     16.0

        The 1.000 at the top is what makes the 16 readable: with no
        modulation the cycle average IS the value.  The mechanism is not
        subtle -- the averaged source injects `4kT <g>` (about half the
        on-state current noise) for the WHOLE period, including the hold
        phase, where the node it injects into is 1 Gohm in parallel with
        100 pF; Hull & Meyer's own condition ("none of the large-signal
        state variables may change significantly over the decay time of
        the impulse response") fails there by six orders (100 ns closed,
        0.1 s open).  So this route is for a mixer's `gm`, a bias-
        dependent shot noise -- not for a switch.  ⚠ A switch's noise IS
        reachable exactly: `covariance` and `oscillator_covariance`
        evaluate `CY` at every step and need no averaging (the switched
        capacitor's held variance reads `kT/C` to 1e-4 there).

        ⚠ THIS IS WHAT `_cy_reduced` REFUSES, DONE INSTEAD OF REFUSED, and
        the literature's answer rather than ours.  Hull & Meyer (1993):
        *"cyclostationary noise sources, such as shot noise, may be modeled
        as MODULATED STATIONARY NOISE SOURCES.  The impulse response that
        is calculated INCLUDES THE EFFECT OF THIS MODULATION.  In the case
        of shot noise, the hypothetical stationary noise source has
        spectral density `S_i = 2q Ibar_c`"* with `Ibar_c` the
        cycle-averaged current.

        So the modulation is carried by the RESPONSE, which `pnoise`
        already computes as `H_l`, rather than by the SOURCES.  Okumura's
        route puts one independent stationary source per timestep interval
        per device -- `p` per device, ~25,000 sources on a real circuit.
        This is ONE per device.  Same physics, `p` times cheaper.

        ⚠ AND ITS CONDITION IS CHECKABLE RATHER THAN A BLANKET REFUSAL:
        *"valid when the impulse response duration is much less than the
        time it takes for the mixer circuit to significantly change its
        state ... NONE OF THE LARGE-SIGNAL STATE VARIABLES MAY CHANGE
        SIGNIFICANTLY OVER THE DECAY TIME OF THE IMPULSE RESPONSE."*

        ⚠⚠ WHICH IS THE OPPOSITE OF HIGH-Q, AND THEY SAY SO: *"high-Q
        filters should be avoided, since they cause the impulse response to
        ring, and thus require a very large value of M."*  So this
        construction degrades exactly where `lambda_2 -> 1` -- the same
        boundary as everything else in this class, arriving from a fourth
        direction.  That makes the two constructions COMPLEMENTARY rather
        than competing: Hull & Meyer for fast-settling circuits, Okumura's
        expensive one for the high-Q case that needs it.  `info` reports
        `|lambda_2|` so the caller can see which regime they are in.

        ⚠ SAMPLED ON THE ORBIT, NOT AT THE OPERATING POINT.  The average
        that matters is over the LARGE-SIGNAL waveform, so `CY` is
        evaluated at every stored state and averaged with the step weights
        -- the same quadrature `diffusion_constant` uses, so the two remain
        comparable.
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

        ⚠⚠ AND ITS SCOPE IS WIDER THAN "AN UNUSUAL CASE" -- IT IS A
        BLANKET REFUSAL OF MOS pnoise.  ⚠ AN EARLIER VERSION OF THIS NOTE
        SAID IT WAS UNREACHABLE BECAUSE `PspMosLongChannel.CY` IS
        IDENTICALLY ZERO.  That was measured on a DEFAULT-CONSTRUCTED
        element: the model has channel thermal and flicker noise, and
        `fnt = 0` by default because "an element built without a card is
        noiseless".  With `fnt = 1` it is nonzero, white, and
        bias-dependent -- so the refusal is REACHABLE from a real device
        today, and `modulated=True` is the route past it.  There is no physically correct MOS noise model whose
        `CY` is state-independent: thermal channel noise is
        `4kT gamma g_d0` with `g_d0` bias-dependent, flicker goes as
        `I_D^AF`, gate shot noise as `2qI_G`, and Mahmutoglu & Demir
        (2015) are explicit that trap rates "depend on the voltages across
        the MOSFET which can considerably vary with time during
        large-signal operation".  So the answer to "will a real device
        pass this check" is already determined, and it is no.

        ⚠ THE ROUTE OUT IS THE CYCLOSTATIONARY CONSTRUCTION, NOT A
        DIFFERENT DEVICE MODEL, and that reorders the roadmap: the
        cyclostationary path is not an enhancement for MOS pnoise, it is
        the PRECONDITION -- for the thermal and shot mechanisms.  ⚠ NOT
        for flicker: Okumura's own construction excludes it (p. 585,
        "cannot be modeled as a cyclostationary process by using this
        method, because it has very long time constants"; verified at the
        source 2026-09-08), by the same long-time-constant physics as the
        trap-rate caveat beside it, and falls back to a stationary flicker
        source.  Hull & Meyer (1993) make it affordable -- one
        stationary source per device at the cycle-averaged current, with
        the modulation carried by the impulse response `H_l` that A1
        already computes -- and their worked example IS shot noise
        modulated by the collector current, i.e. exactly this case.  Their
        condition is checkable rather than a blanket refusal, and it fails
        in the familiar direction: a ringing impulse response breaks it,
        so it degrades as `lambda_2 -> 1`.

        ⚠ SECOND-ORDER CONSEQUENCE, WORTH KNOWING BEFORE THE MODEL LANDS.
        This same check is what keeps the Ito/Stratonovich choice out of
        reach (`CY = GG^T`, so a state-dependent `CY` is a state-dependent
        `G`).  Relaxing it for MOS makes the two interpretations diverge,
        and Demir's escape -- "the noise signals are small compared with
        the deterministic signals" -- may NOT carry for trap noise: a trap
        occupancy is a two-state Markov chain rather than a small
        perturbation of a large signal, and the same paper says the state
        dependence "in fact makes the equation nonlinear".  The tell would
        be a discrepancy in a MEAN but not in a variance.
        """
        irn = pss.irefnode
        fp = pss.factored_period()
        ## ⚠ THREE STATES ON THE ORBIT -- the first used to be the ZERO
        ## VECTOR, which is on the orbit only by accident, and a linear
        ## time-invariant RC held by a DC clock was refused as
        ## cyclostationary because a switch model read `goff` at v(ck) = 0
        ## (found by an external reference-simulator cross-check, 2026-09-05).  The third
        ## probe is now the stored state half a period in.
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
    ## `f0`.
    ##
    ## ⚠ TIGHTENED BY A7.  This used to be the conditioning floor being
    ## accepted, because `sigma_min` of the plain operator falls LINEARLY
    ## with the distance and everything nearer was unusable.  The deflated
    ## solve removes that: its conditioning is FLAT (measured 2.04e-01 from
    ## 0.3 down to 1e-9 of `f0`), so the only remaining reason to refuse is
    ## the physical one -- at an EXACT harmonic `1/(1 - alpha)` is a
    ## division by zero and the response is genuinely unbounded.  So the
    ## guard now excludes only what has no finite answer, not what was
    ## merely hard to compute.
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
        circuit's unit multiplier makes singular.  MEASURED on van der Pol,
        `sigma_min(I - exp(-jwT) M)`:

            offset/f0    0     0.25    0.5    0.75    1      2      3
            sigma_min   2.8e-11 0.51   0.65   0.51  2.8e-11 2.8e-11 2.8e-11

        and LINEAR in the distance to the nearest one -- 2.5e-1, 2.6e-2,
        2.6e-3, 2.6e-4 at 0.9, 0.99, 0.999, 0.9999 of the way there.  The
        same sweep on a DRIVEN ladder never drops below 0.78: no unit
        multiplier, no singularity, harmonics included.

        ⚠ AND IT IS PHYSICS, NOT CONDITIONING.  A perturbation at a
        harmonic is a perturbation along the orbit, and an oscillator's
        response to that is unbounded phase drift -- there is no bounded
        periodic answer to return.  So this refuses rather than tightening
        a tolerance, and says which quantity to ask for instead.
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
        """GMRES, judged by its RESIDUAL rather than by its status flag.

        ⚠ SCIPY REPORTS BREAKDOWN ON SYSTEMS IT HAS ALREADY SOLVED.  These
        operators are `2m x 2m` and often tiny, so the Krylov space is
        exhausted in a handful of steps; the next vector is then numerically
        zero, which is a LUCKY breakdown -- the solution is exact -- and it
        comes back as `info = 4` all the same.  Trusting the flag turns an
        exact answer into a `RuntimeError`, which is what it did for AM/PM
        at small offsets.

        So the residual decides.  A genuine failure still fails, and it
        fails with the residual quoted, because the real cause near a
        harmonic is that the operator is nearly singular there and no
        tolerance will fix it.
        """
        ## ⚠ NOW OUR OWN ARNOLDI-GMRES, which returns the residual as its
        ## verdict instead of a status flag that has to be overruled.
        ## The workaround below survives as the TOLERANCE decision; what
        ## has gone is the second opinion about whether the solve failed.
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

        ⚠ THE TRAP IS THAT NOTHING ELSE WOULD OBJECT. The Lyapunov
        recursion, `diffusion_constant` and eq (22)'s collapse all read
        `CY` at ONE frequency and treat it as the noise intensity at every
        frequency; a coloured source folded that way returns a plausible
        number, not an error (A4d names exactly this shape). Detected by
        evaluating the reduced `CY` at two frequencies -- colour is
        frequency dependence, bias dependence is what `_cy_reduced`
        refuses separately.
        """
        w1 = 2.0 * np.pi / float(pss.period)
        ## ⚠ ONE state, two frequencies: the colour question is separable
        ## from the bias question, and asking it through `_cy_reduced`
        ## refused every MODULATED source before the covariance routes
        ## (which evaluate `CY` per step and handle modulation exactly)
        ## could reach it.
        ## ⚠⚠ PER ENTRY, AND AT MORE THAN ONE STATE (2026-09-15; reported by a
        ## peer session, reproduced).  This compared the difference against
        ## `1e-9 * max|CY|` over the WHOLE matrix, at `x_last` only.  Once the
        ## PSP gate resistor carried its white 4kT/rg (`d0e7e4e`), a low-rg
        ## device put `1.27e-20` in `CY` and a drain's flicker colour
        ## (`2.5e-30` on its own `1.3e-27`) passed that global threshold --
        ## `covariance` folded 1/f at w0 on a sample-and-hold and said
        ## nothing.  And a source white at the period boundary (a switch OFF
        ## at t = 0) but coloured elsewhere passed at any scale.  So each
        ## entry is judged against ITS OWN magnitude, at `x_last` and at
        ## states spread over the orbit.  Exact zeros and white entries give
        ## identical values at both frequencies, so neither can fire.
        ## `_cy_at` takes a REDUCED state: the reference row is removed here.
        _xl = np.asarray(pss.factored_period().x_last, dtype=float).ravel()
        _states = [_xl]
        _wf = getattr(pss, 'waveform', None)
        if _wf is not None:
            _W = np.delete(np.asarray(_wf[1], dtype=float), pss.irefnode,
                           axis=0)
            for _k in sorted(set(np.linspace(0, _W.shape[1] - 1,
                                             8).astype(int))):
                _states.append(_W[:, _k])
        coloured = False
        for _xr in _states:
            c1 = self._cy_at(pss, w1, _xr)
            c2 = self._cy_at(pss, 10.0 * w1, _xr)
            den = np.maximum(np.abs(c1), np.abs(c2))
            if np.any(np.abs(c1 - c2) > 1e-9 * den):
                coloured = True
                break
        if coloured:
            raise NotImplementedError(
                'PAC.%s: a noise source in this circuit is COLOURED (its CY '
                'differs between w0 and 10 w0), and this routine assumes '
                'white sources -- it would fold CY at one frequency as if '
                'it held at every frequency and return a plausible wrong '
                'number. Use the frequency-resolved surfaces (pnoise, '
                'phase_psd/coloured_diffusion), or the white-through-filter '
                'form of the source.' % what)

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
        one implementation, so the pair cannot drift apart in the way that
        `diffusion_constant` and `covariance` once did over exactly this
        factor of two.
        """
        self._refuse_coloured(pss, what)
        fp = pss.factored_period()
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
        ## period.  The Lyapunov accumulation was already per step; hoisting
        ## a single `CY` out of it put a MODULATED source (a switch's
        ## `4kT g(t)`, a MOS channel's `4kT gamma gd0(t)`) outside the
        ## formulation rather than outside the accuracy, and the
        ## cyclostationarity refusal in `_cy_reduced` then closed the door
        ## on exactly the circuits whose noise is the point (found by the
        ## reference-simulator cross-check, 2026-09-05).  Evaluated at the state
        ## the step's companion was factored at (the implicit step's own
        ## solution); the colour refusal still applies -- colour is a
        ## different axis.
        w0 = 2.0 * np.pi / float(fp.T)
        _W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                       pss.irefnode, axis=0)
        cys = [np.real(self._cy_at(pss, w0,
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
            ## zero-stability bound on an event grid, 2026-09-21) has no
            ## third alpha: the pair map holds with it zero -- see
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

        The period map on that pair is the plain product of the `A_k`
        with NO re-seeding of `iq` at the boundary.  ⚠ THAT IS A DIFFERENT
        OBJECT FROM `fp.matvec`, deliberately: the shooting SOLVE re-seeds
        the companion at each period start (the manufactured opener, B16),
        but the discretised noisy system does not, and the covariance is
        a property of the latter.  The tie between the two is exact and is
        the gate: the pair product applied to `(x, 0)` and read out on `x`
        IS `fp.matvec`.

        ⚠ `oscillator_covariance` IS REFUSED FOR TRAP-PLAIN, with the
        reason: it borders `I - M kron M` with `ppv()`'s null vectors,
        which are width `m` on the plain path, and the pair map is
        `2m x 2m`.  The pair's own null vectors would be needed, with a
        normalisation this record has been burned on twice today
        (`floquet_modes`' state-block scale, `ppv()`'s `v . xdot`).
        Named rather than approximated; euler-plain and gear both work.
        """
        m = pss.cir.n - 1
        hs = np.diff(np.asarray(fp.times, dtype=float))
        ## ⚠ `CY` PER STEP, AT THE STEP'S OWN STATE -- not one `CY` for the
        ## period.  The Lyapunov accumulation was already per step; hoisting
        ## a single `CY` out of it put a MODULATED source (a switch's
        ## `4kT g(t)`, a MOS channel's `4kT gamma gd0(t)`) outside the
        ## formulation rather than outside the accuracy, and the
        ## cyclostationarity refusal in `_cy_reduced` then closed the door
        ## on exactly the circuits whose noise is the point (found by the
        ## reference-simulator cross-check, 2026-09-05).  Evaluated at the state
        ## the step's companion was factored at (the implicit step's own
        ## solution); the colour refusal still applies -- colour is a
        ## different axis.
        w0 = 2.0 * np.pi / float(fp.T)
        _W = np.delete(np.asarray(pss.waveform[1], dtype=float),
                       pss.irefnode, axis=0)
        cys = [np.real(self._cy_at(pss, w0,
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
        if pair and what == 'oscillator_covariance':
            raise NotImplementedError(
                'PAC.oscillator_covariance: the trapezoidal plain path\'s '
                'per-step state is the pair (x, iq), so its period map is '
                '2m x 2m, and the bordered solve needs THAT map\'s null '
                "vectors -- ppv()'s are width m. Not built (see "
                '_lyapunov_pieces_plain). Use method=\'euler\' for a plain '
                "width-m reference, or method='gear'.")
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
            ## per step, so the un-reset pair carries a marginal mode --
            ## the `(-1)^n` obstruction this file records for every
            ## formulation that keeps `iq` across a period (measured here:
            ## LinAlgError on a driven RLC).  The shooting solve is
            ## well-posed because the manufactured opener re-seeds `iq` at
            ## zero; the covariance's period map must do the same.  With
            ## `iq` zeroed at the start, the x->x block of the product IS
            ## `fp.matvec` (tied to 1e-12 above), and the map on the pair
            ## is the product applied to `(x, 0)`.
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
        (Roemisch & Winkler; confirmed by the naive two-stage scheme coming
        out 27% biased on kT/C).  Van Loan evaluates it exactly: the
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

        Verified against kT/C at second order on R||C (ODE) and VS-R-C
        (DAE), matching an independent measurement to the digit; the
        stationary error is the METHOD's O(h^2), NOT machine zero (a
        machine-zero kT/C would mean a method-consistent `Q = P(1-A^2)`
        fudge that corrupts the transient covariance).
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

    def _lyapunov_pieces_stage(self, pss, fp, what):
        """`_lyapunov_pieces` for a Runge-Kutta stage method's Floquet source
        (Radau IIA, TR-BDF2, ESDIRK).

        Identical in structure to `_lyapunov_pieces_trbdf2`: the per-step
        transition `A_n` is the stage step map (dense, `m x m`, via
        `_monodromy_matvec_stage` one step at a time) and the per-step
        injection `Q_n` is the stage injection (`_stage_injection`) -- the
        source enters every STAGE; the end-of-step DAE-projected VAN LOAN
        integral (`_vanloan_step_injection`) read a switch's held variance
        0.876 kT/C at 400 points (O(h)) and stays the fallback.  State width
        `m`, so `n = m`.

        ⚠ THE VAN LOAN INJECTION IS EXACT, THE METHOD SETS ONLY THE
        PROPAGATION.  Under it the covariance still converges to the
        stationary target (kT/C on an RC) at the injection's O(h^2), not at
        Radau's O(h^5): Van Loan already integrates the step exactly, so
        refining the grid gains on the recursion's discretisation of a
        continuous Lyapunov flow, which the higher-order transition does not
        change.
        """
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
            CYn = self._cy_at(pss, w0, xk)
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
        the source entering EVERY stage (2026-09-15), or None when the
        tableau has a non-positive weight (then the caller keeps Van Loan).

            Q_k = sum_i T_i (CY(Y_i) / (2 h b_i)) T_i^T,
            T_i = d x_{k+1} / d u_i   through the method's own stage solve

        with `CY` at the STAGE states `Y_i`.  White noise over the step is
        the method's quadrature `h sum_i b_i u_i` with independent stage
        samples of variance `CY/(2 h b_i)`, so the increment's variance is
        `h CY/2` -- the diffusion -- and each sample reaches the step's end
        through the stage equations exactly as a stage source does.

        ⚠⚠ WHY (measured on a switched capacitor, Ron 1 k, 100 pF, 100 kHz).
        The Van Loan injection freezes `C`, `G`, `CY` at the END of the step,
        so across the switch-off edge it integrates the injection with the
        OFF conductance: held variance 1 - 0.876 / 0.934 / 0.966 kT/C at
        400 / 800 / 1600 points under radau (first order), 0.87 under
        trbdf2.  With the stage injection radau reads 1.3e-7 off at 400
        points and trbdf2 2.5e-3 (second order, 4.1x per doubling); a
        b-weighted Van Loan at the stage states read 5.5e-4 / 4.6e-3.  On a
        constant-operating-point RC it is not an exactness fit: the error is
        nonzero and falls with the grid (radau 1.4e-9 / 4.4e-11 / 1.4e-12,
        trbdf2 2.6e-4 / 6.3e-5 / 1.6e-5 at 100 / 200 / 400 points).
        ⚠ Needs stiff accuracy (`x_{k+1} = Y_s`) and positive weights:
        radau and trbdf2 qualify.  ⚠⚠ A TABLEAU WITH A NON-POSITIVE WEIGHT
        (ESDIRK43: `b = 0.158, 0, 0.187, 0.681, -0.275, 0.25`) takes the Van
        Loan injection at the STAGE states instead, averaged with POSITIVE
        trapezoid weights over the stage abscissae in time.  Measured on the
        same sampler, 1 - held at 100 / 200 / 400 points: end-of-step Van Loan
        0.367 / 0.224 / 0.124; this 2.1e-2 / 3.4e-3 / 3.7e-4 (tracking and a
        constant-operating-point RC identical to the end-of-step form -- the
        weights sum to one); equal-variance stage samples, the other
        tableau-independent candidate, 5.3e-2 / 3.6e-2 / 2.0e-2 and 6.8e-2 off
        on the RC -- rejected.
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
                CYi = np.real(np.asarray(self._cy_at(pss, w, yi), dtype=complex))
                Q += wts[i] * self._vanloan_step_injection(
                    np.asarray(pss._C_at(yi), dtype=float),
                    np.asarray(pss._G_at(yi), dtype=float), CYi, h)
            return 0.5 * (Q + Q.T)
        Q = np.zeros((m, m))
        for i in range(s):
            yi = np.delete(np.asarray(states[k * s + i], dtype=float), irn)
            CYi = np.real(np.asarray(self._cy_at(pss, w, yi), dtype=complex))
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
        DAE'S OWN DERIVATIVE at the node's state (2026-09-22): on the
        differential rows ``C(x) xdot = -(i(x) + u(t))``, on the algebraic
        rows (a zero row of `C`) the differentiated constraint ``G(x) xdot
        = -du/dt``; one small solve per node, exact for the discrete state
        and independent of the step.  The rate converts a node's motion in
        time into a state change: see `_fixed_time_event_columns`.

        ⚠ IT WAS A THREE-NODE PARABOLA (one-sided at a landed event), and
        that cost 5 % on the comparator oscillator's collapsed `c` node
        inside its 10 ns ON phase, where a 7 ns step cannot fit a parabola
        to a 10 ns exponential -- the one node of the staged-oscillator
        sideband test that had to be excluded.  Against the exact
        piecewise-linear rate the DAE form reads 1e-15 at every node
        outside the windows (the stencil 10.6 % at worst, 0.45 % at the
        node that had to be excluded); on the driven jitter sampler the
        sawtooth's rate is its slope exactly and the held capacitor's its
        leak.  The stencil is kept only as the fallback where the
        assembled matrix is singular (an index above one), and says so."""
        ts = np.asarray(pss.waveform[0], dtype=float)
        X = np.delete(np.asarray(pss.waveform[1], dtype=float),
                      pss.irefnode, axis=0)
        N = len(ts) - 1
        m = X.shape[0]
        out = np.zeros((N + 1, m))
        analysis = getattr(pss.par, 'analysis', None)
        ok = True
        for j in range(N + 1):
            ## ⚠ NODE 0 IS EVALUATED AS NODE N (2026-09-23).  A source's
            ## derivative at exactly its start (`VSin` clamps `t - td` at 0:
            ## SPICE's rule for a transient) is the LEFT one -- 0 -- where the
            ## periodic steady state, t = 0 == t = T, has the right one: the
            ## rate read 0 for an exact 6283 V/s at node 0 on a driven RC.
            ## Harmless to the consumers (the node's own motion `tau_0` is 0)
            ## and wrong as a function; node N is the same point.
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
        """The three-node parabola `_orbit_rate` used until 2026-09-22 --
        one-sided AT a landed event and at the node after one, central
        elsewhere, periodic at the ends.  Kept as the fallback."""
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
        """The BORDERED Lyapunov closure on a staged solve (2026-09-22,
        events phase B), or `None` when the solve is not staged.

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
        a gain the three collocation points invent -- a comparator-jitter
        sampler (a noisy threshold, a hold capacitor on a second ramp)
        read 3.8x its analytic held variance `(s_2/s_1)^2 kT/C_n`.  With
        the crossing conditions pinned at both window edges the bordered
        system cancels that internal sensitivity, which is why the
        moving events must be unknowns of the noise problem too.

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
        its gear twin, which has no columns -- warned, unbordered)."""
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
        ## gear's PAIR form (item 2 of the list after E7, 2026-09-22): the
        ## state is (x_j, x_{j-1}), the event row acts on the first block,
        ## the per-node column of node j is the pair (Pk_j, Pk_{j-1}) and
        ## the map to node j the pair of `P_nodes` rows; the samples come
        ## out as pair covariances, as the plain gear path returns them
        W = [np.pad(np.asarray(w, dtype=float).ravel(), (0, n - m)) for w in ev['W']]
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
        d = np.zeros((N, K, n))
        for k, nd in enumerate(nodes):
            r = W[k].copy()
            for j in range(nd - 1, -1, -1):
                d[j, k] = r
                r = r @ As[j]
        Z = np.zeros((N + 1, n, K))
        Kf = np.zeros((N + 1, n, n))
        D = np.zeros((K, K))
        for j in range(N):
            Z[j + 1] = As[j] @ Z[j] + Qs[j] @ d[j].T
            Kf[j + 1] = As[j] @ Kf[j] @ As[j].T + Qs[j]
            D = D + d[j] @ Qs[j] @ d[j].T
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

    def event_jitter(self, pss):
        """The noise-driven JITTER of every landed crossing of a staged,
        driven solve (2026-09-22): ``sigma`` in seconds per crossing, and
        the crossings' covariance in fractions of the period.

        The bordered Lyapunov closure (`_event_closure`) already carries
        it: the crossings move as ``dtheta = (dtheta/dx_0) dx_0 - Gt^-1
        sum_j d_j w_j`` -- the stationary state at the period start
        (covariance `K_0`, from the previous periods' noise) and this
        period's per-step injections, independent of each other -- so
        ``Cov(dtheta) = dth K_0 dth^T + Gt^-1 D Gt^-T`` with ``D = sum_j
        d_j Q_j d_j^T``.  Measured on the comparator-jitter sampler
        (`_jitter_sampler`: a sawtooth of slope s_1 crossing a threshold
        node with kT/C_n of noise): the turn-off crossing's sigma is
        ``sqrt(kT/C_n) / s_1`` -- 11.5873 ps measured against 11.5844 ps
        analytic, 1.0002, and flat at 100 / 200 / 400 points -- the same
        crossing motion that gives the held capacitor its
        `(s_2/s_1)^2 kT/C_n`.  The reset edges of that fixture read
        0.6437 ps, the threshold node's own faster slope there.

        Returns ``{'sigma': (K,) s, 'cov_fraction': (K, K), 'fractions':
        (K,) the crossings' positions, 'nodes': (K,) their grid nodes}``.
        An oscillator's crossings diffuse without bound with its phase;
        that is `oscillator_covariance`'s object, and this refuses one.
        Every source of the circuit is in it together; a per-source
        split is the per-source `Q_j`, which `_lyapunov_pieces` does not
        keep."""
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
        cov = 0.5 * (cov + cov.T)
        T = float(pss.period)
        return {'sigma': np.sqrt(np.clip(np.diag(cov), 0.0, None)) * T,
                'cov_fraction': cov,
                'fractions': np.asarray(pss._state_event_fracs, dtype=float).copy(),
                'nodes': list(pieces['nodes'])}

    def covariance(self, pss, samples=False):
        """The periodic (cyclostationary) state covariance — DRIVEN circuits.

        ⚠ A GRID CHOSEN FOR `kT/C` IS NOT A GRID FOR THE PROFILE.  With
        `CY` per step (2026-09-05) a switched capacitor's HELD variance
        reads `kT/C` to 1e-4 at 1600 points and converges at better than
        second order, while the TRACKING phase sits at this routine's
        O(h/tau) floor -- 4% out at 800 points where the held value is
        already 1.6e-4 -- and both agree with a reference simulator's sampled pnoise at
        matched instants to 1e-3 (0.99878 track, 0.99915 edge, 0.99999
        hold).  ⚠⚠ **AN EARLIER READING OF THAT TRACKING NUMBER IS
        WITHDRAWN (2026-09-18).**  It said: "the tracked variance is 0.957
        kT/C, NOT kT/C -- a sinusoidal clock holds the switch at full `gon`
        only instantaneously, so the capacitor is never in equilibrium with
        `Ron`; both tools agree on that independently."  That is a
        DISCRETISATION FLOOR described as physics.  Within this fixture's
        own definition the switch contributes `g V` and `white_noise(4 kT
        g)` with the SAME `g`, so the variance obeys `dV/dt = -2(g/C)V +
        2kT g/C^2`, for which `V(t) = kT/C` is an EXACT solution at every
        instant and for ANY `g(t)` -- fluctuation-dissipation, and the
        periodic solution is unique because `g > 0` contracts.  The
        capacitor IS in equilibrium throughout.  Measured against that
        exact profile (a stiff solve to periodic steady state, closing
        error 0.0): the tracking value reads 0.7398 / 0.8467 / 0.9157 /
        0.9556 kT/C at 200 / 400 / 800 / 1600 points -- errors 0.260,
        0.153, 0.084, 0.044, HALVING per doubling, i.e. first order,
        converging to kT/C.  ⚠ And the cross-tool agreement does not
        rescue the claim: two tools discretising the same period at
        comparable step counts agree about a SHARED artefact, which is
        exactly why agreement between implementations cannot establish a
        LIMIT -- the same lesson the edge-jitter work paid for when
        cross-family agreement could not show a number was right.
        ✅ CONSISTENCY IS SEPARATELY CONFIRMED, and it is what this routine
        should be judged on: against a closed-form time-varying reference
        built by breaking the fluctuation-dissipation balance (an extra
        white source not tied to `g`, giving a profile that spans 15x over
        the period), `covariance` converges to the exact continuous answer
        at FIRST order in both phases -- hold 1.07e-2 -> 1.09e-3 and track
        2.60e-1 -> 4.44e-2 over 200..1600 points, ratios 2.07 and 1.90.
        That is the O(h) its piecewise-constant injection predicts, and it
        converges to the RIGHT limit.  ⚠ The "held converges at better
        than second order" above is a property of the kT/C fixture, where
        the exact profile is a CONSTANT and the leading terms cancel: with
        the balance broken the held value converges at first order too.
        ⚠⚠ **AND THE kT/C ARGUMENT IS SPECIFIC TO A FIXTURE WHOSE NOISE IS
        TIED TO ITS OWN CONDUCTANCE** -- the caveat is the reference suite's
        and it is right.  It is a theorem about THIS circuit only because
        the element's own definition contributes `g V` and `white_noise(4 kT
        g)` with the SAME `g`.  A REAL DEVICE WOULD BREAK IT -- on a PSP
        switch the measured `sid/(4kT g)` runs 1.09 during conduction to
        3.17 through turn-off, so fluctuation-dissipation would not balance
        and that circuit's tracking limit need not be kT/C.  ⚠ But that is
        a caveat about OTHER fixtures, not about this one: the reference
        side runs the SAME behavioural switch (`pcswitch.va`, the same `g`
        and the same `white_noise(4 kT g)`), which is the point of the
        fixture, so the theorem applies to both sides and there is no real
        device near it.
        ⚠⚠ A REFERENCE-SIDE READING OF ~0.974 kT/C FOR THE TRACKING LIMIT
        IS WITHDRAWN (2026-09-18, by the side that made it).  It came from
        an Aitken extrapolation of 0.95489 / 0.96723 / 0.97128 / 0.97239,
        whose error ratios DECELERATE -- 1.38, 1.14, 1.04.  A ratio heading
        to ONE is a sequence that has stopped moving, and Aitken on a
        stalled sequence returns approximately where it stalled rather than
        a limit; a converging first-order ladder heads to TWO, as the one
        above does (1.70, 1.82, 1.90).  So that number is a floor in the
        reference's own sampled-noise integration, and the limit for this
        fixture remains kT/C by the argument above.

        Returns `K0`, the covariance at `t = 0`; with `samples=True`,
        `(K0, [K_j])`, the covariance at every step, which is the
        time-varying statistic this exists to produce.

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
        `I - M` — measured 3.1e-11 against 3.8e-11 — and the covariance
        does not settle, it GROWS.  Variance linear in `t` is a random
        walk, which is phase diffusion, which is the linewidth.  Demir 2002
        gives the physical counterpart: an oscillator's output noise is
        STATIONARY, not cyclostationary, because "noisy autonomous systems
        cannot provide a perfect time reference".  There is no
        cyclostationary object there to compute, and `oscillator_spectrum`
        is the right route instead.

        ⚠ `CY/2` IS THE ONE-SIDED-TO-TWO-SIDED CONVERSION AND IT IS NOT
        COSMETIC.  `CY` is a one-sided density (a resistor's `4kT/R`), so
        the per-step injection is `Q_j = Jf_j^-1 (CY_j / 2h_j) Jf_j^-T`.
        MEASURED against `kT/C` — exact, famously independent of `R` — on
        an RC circuit, with and without the half:

            npts        100      200      400      800
            CY          1.861    1.928    1.963    1.981
            CY/2        0.931    0.964    0.982    0.991

        The full-`CY` column converges to 2 and the halved one to 1, so the
        factor is settled by the measurement rather than by argument.  The
        residual halves per grid doubling — O(h), first order, which a
        piecewise-constant approximation to white noise is.

        ⚠ AND THE GRID MUST RESOLVE THE NOISE BANDWIDTH, which is a real
        precondition rather than an accuracy note.  The first attempt at
        that gate read 0.517 because the RC pole at 159 kHz sat ABOVE the
        grid's 100 kHz Nyquist: the discrete system genuinely does not
        carry the noise the continuous one does.  A `kT/C` that comes back
        low is the grid, not the code.

        ⚠ COST: the solve has `(2m)^2` unknowns and is dense here, so it is
        `O(m^4)`.  Small circuits only until that is replaced.
        
        ⚠ ON A STAGED SOLVE (`state_events=True`) THE CLOSURE IS BORDERED
        (2026-09-22): the noise moves the landed crossings, and the plain
        closure on such a solve is wrong by O(1), not merely incomplete --
        see `_event_closure` for the algebra and the comparator-jitter
        sampler that measured 3.8x its analytic held variance unbordered
        and 0.999 bordered.  Samples are the covariance at FIXED times.
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
        As, Qs, K1, M, m, n = self._lyapunov_pieces(pss, 'covariance')
        ## a staged solve closes on the TOTAL monodromy with the events'
        ## noise-driven motion in the injection -- see `_event_closure`
        bordered = self._event_closure(pss, As, Qs, M, m, n)
        if bordered is not None:
            M, K1, _samples, _pieces = bordered
        S = np.eye(n * n) - np.kron(M, M)
        K0 = np.linalg.solve(S, K1.reshape(-1)).reshape(n, n)
        K0 = 0.5 * (K0 + K0.T)
        if not samples:
            return K0
        if bordered is not None:
            return K0, _samples(K0)
        seq, K = [K0], K0
        for A, Q in zip(As, Qs):
            K = A @ K @ A.T + Q
            seq.append(0.5 * (K + K.T))
        return K0, seq

    def sampled_noise(self, pss, output, times, freqs, maxsidebands=None,
                      tail=False):
        """The one-sided PSD of the SAMPLE SERIES `y(t0 + kT)` -- DRIVEN
        circuits, white AND coloured sources.  2026-09-15.

        Returns `S` of shape `(len(times), len(freqs))` in `output`'s units
        squared per Hz, for `0 < f <= f0/2`.  `sum` over the band of `S` is
        the variance at `t0` that a sampler sees; see `sampled_variance`.
        The instants actually used (the nearest period-grid points) are left
        in `self.sampled_instants`.  ⚠ When comparing at a round instant,
        pass the grid's own times (`pss.factored_period().times`): the grid
        need not have the step count the `timestep` suggests (T/1000 can give
        999 steps), and with a 1/f source the held variance moves measurably
        between neighbouring grid points (peer report).

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
        A joint root of the summed `CY` made independent sources
        non-additive (measured +9.3 % of the held variance, white + flicker
        in one element).  ⚠ Mahmutoglu & Demir (TCAS-I 62(4), 2015) show
        that a SWITCHED MOSFET's trap (1/f) noise is whitened below the
        switching frequency and that a modulated-stationary 1/f model
        over-predicts it; the physical fix needs trap states in the device
        model, which is outside this analysis.  An agreement with another
        tool on this quantity is agreement on the convention.  A flicker
        spectrum is also singular at DC: keep `f` away from 0 (the band
        integral takes an explicit `fmin`).

        MEASURED (switched capacitor, Ron 1 k, Roff 1 G, 100 pF, 100 kHz,
        gear, 400 points): the white held variance equals `covariance`'s
        at the same instant to 1e-6 (0.998721 kT/C) and the tracking value
        to 1.6e-4 (0.8466, covariance's own O(h/tau) floor); the seeded row
        at `t0 = 0` equals `adjoint_transfer_row` to 1.5e-16; on an LTI
        circuit the series sum equals the fold of `pnoise` over the same
        sidebands to 1e-10, white and 1/f.  Each gate was checked to fail:
        the source samples one step early read 1.30 kT/C, sidebands cut to
        N/8 read 0.77 in track.

        ⚠ STAGE METHODS (radau, trbdf2) RUN NATIVELY, 2026-09-15: the source
        enters every stage, so the sensitivities are read at the stage
        abscissae and `CY` at the STAGE states (one re-traversal of the
        orbit, cached).  Measured on the sampler: the LTI series sum equals
        the pnoise fold to 1e-10 under both; the held variance reads
        1.000000 (radau) / 0.999979 (trbdf2) kT/C at 400 points and stays
        there at 800, the tracking value converges at first order (radau
        0.949 -> 0.975, trbdf2 0.916 -> 0.957).  ⚠ `covariance` under these
        methods WAS no reference at a switching edge until its injection
        moved to the stages too (`_stage_injection`, 2026-09-15): the Van
        Loan injection frozen at the END of each step read the held variance
        0.876 / 0.934 (radau, 400 / 800 points); it now reads 1.3e-7 off at
        400.  GLM period maps are refused.

        ⚠ TIME AVERAGE, PER FREQUENCY: the mean over `t0` of this PSD is
        the fold of the time-averaged PSD, `sum_k pnoise(|f + k f0|)`, over
        EVERY output band to the grid's Nyquist -- measured 4.4e-5 at 100
        points (a fold cut at |k| <= 10 left 0.4 %).

        ⚠ TWO THINGS LIMIT A HELD VARIANCE, AND THEY ARE NOT THE SAME THING
        (2026-09-21, from a peer's kT/C ladder on a half-on switch).  The
        fold stops at `|n| <= maxsidebands`, so the source spectrum beyond
        `F = (L + 1/2) f0` is dropped: for a Lorentzian that is the tail
        `1 - (2/pi) atan(F/fc)` ~ `(2/pi) fc/F`, first order in the sideband
        count.  `tail=True` adds, per instant and series frequency, the
        `1/nu^2` extrapolation of the two OUTERMOST covered sidebands over
        the uncovered half-lines -- `A = nu_L^2 dens_L`, `A / (f0 (F +
        f0/2))` each side -- which is exact for a single pole and correct to
        `(fc/F)^2` in general; nothing is fitted.  AND the covered sidebands
        are computed on the grid: at `omega h ~ 3` per step (the top bands
        of a fold to the grid's Nyquist) gear's discrete transfer is far
        from the continuous one.  Measured on the LTI limit of the switched
        capacitor (g = 0.5 mS, C = 100 pF, fc = 796 kHz, f0 = 100 kHz), the
        deficit against the pure tail at a FIXED 100 sidebands: 3.03 /
        1.93 / 1.33 / 1.12 / 1.05 at 204 / 400 / 800 / 1600 / 3200 points --
        the peer's "coefficient 3.0" was this, at a fixed M/npts, not a
        property of the kernel.  A run whose top sideband sits above
        `omega h = 1` is WARNED (`SAMPLED_RESOLUTION_WARN`): the remedy is
        more points; `tail=True` then closes what lies beyond the covered
        edge, which needs that edge well above the spectrum's corner (edge
        at 12 fc: radau 0.948 -> 0.998 x kT/C at 204 points, gear 0.947 ->
        0.995 at 3200; edge at 0.8 fc: 0.71, the 1/nu^2 form is wrong
        inside the corner).
        """
        return self._sampled_series(pss, output, times, freqs, maxsidebands,
                                    tail=tail)

    def sampled_variance(self, pss, output, times, fmin, fmax,
                         points_per_decade=40, maxsidebands=None, tail=False):
        """The variance at sampling instants over the SERIES band
        `[fmin, fmax]`, `0 < fmin < fmax <= f0/2`: the trapezoidal integral
        of `sampled_noise` on a log grid of `points_per_decade`, nothing
        added below `fmin`.  Returns an array over `times`.

        ⚠ `fmin` AND `fmax` ARE REQUIRED.  With a 1/f source the integral
        grows as `ln(fmax/fmin)` (measured 0.0010 kT/C per decade, flat from
        1e-4 to 1e-1 f0, on a switched capacitor with a constant 1/f source)
        and has no limit at `fmin -> 0`; with white sources only, the band
        removes `fmin/(f0/2)` of the full variance because the series PSD is
        flat.  Nothing is extrapolated into `[0, fmin]`.
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
        from scipy.integrate import trapezoid
        return trapezoid(S, fs, axis=1)

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

        MEASURED 2026-09-16 on a one-stage linear fixture built so the answer
        is known exactly (`tau = RC = T`, LTI noise path, so `rho_k = e^-k`):
        the transform reproduces `exp(-k)` at **1.0005 at every lag** k = 1..4
        (0.368067 / 0.135404 / 0.049812 / 0.018325 against 0.367879 /
        0.135335 / 0.049787 / 0.018316).  A Monte Carlo over 1176 noisy
        crossings -- no PSS, no adjoint, no spectrum anywhere in it -- agrees
        within 1 sigma at every lag it can resolve (0.346404 / 0.131894 at
        k = 1, 2, i.e. 0.74 and 0.12 sigma), `sigma_t` matching at 0.9901.

        ⚠ `dc_rectangle` EXTRAPOLATES, WHICH THE REST OF THIS FAMILY REFUSES
        TO DO.  `S` is known only on `[fmin, fmax]`; adding `S(fmin)*fmin`
        assumes the series PSD is FLAT below `fmin`.  That is exact for white
        sources and WRONG for `1/f`, where the integral has no limit as
        `fmin -> 0` (see `sampled_variance`).  Hence off by default.
        MEASURED: with it, `rho_k` is unchanged over a 100x range of `fmin`
        (50 -> 0.5 Hz, identical to four digits); without it, `rho_k` drifts
        with `fmin` exactly as truncation should and the drift GROWS with `k`
        (0.9889 of analytic at k = 4, fmin = f0/2e4, recovering to 0.9994 at
        f0/2e5).  ⚠ With a coloured source, LOWER `fmin` -- do not reach for
        the rectangle.

        ⚠ `slew` IS A FINITE DIFFERENCE ON THE PSS GRID, central about the
        instant actually used, and it converges at FIRST order (measured
        0.9924 / 0.9962 / 0.9981 of the analytic slope at 200 / 400 / 800
        points).  It is returned so a caller can check it rather than trust
        it; every metric here is inversely proportional to it.

        ⚠ THE INSTANT IS THE CALLER'S, deliberately.  This does not hunt for a
        threshold crossing: a threshold inferred from a simulated record can
        be biased by startup, which moves the crossing off the steepest point
        -- that cost half of an apparent deficit before it was caught (A8).
        Pass the instant you mean; `instant` in the result is the grid point
        actually used.

        Returns a dict: `sigma_t`, `rho` (k = 1..kmax), `k_cycle`
        (k = 1..kmax), `cycle_to_cycle`, `slew`, `R` (k = 0..kmax), `instant`.
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
        ## ⚠⚠ TWO EARLIER GUARDS HERE WERE WRONG, THE SECOND BECAUSE OF A
        ## MISREAD NUMBER.  `slew == 0.0` is unreachable on a grid, so it was
        ## no guard at all.  Replacing it with a threshold RELATIVE to the
        ## steepest slope looked right only because I had read the peak's
        ## central DIFFERENCE (-8.86e-07) as a SLOPE: over 2h = 5e-9 that is
        ## -1.767e+02, which is 3.58e-04 of the steepest 4.94e+05 -- not the
        ## 1.8e-12 I inferred.  And on a 400-point grid 3.58e-04 is the
        ## SMALLEST ratio any instant can have (a grid point never lands
        ## exactly on the peak), so a relative threshold only ever describes
        ## the grid and moves with N.
        ##
        ## What does not move with the grid is the LINEARISATION this whole
        ## family rests on: a crossing is displaced by delta_y/slew only while
        ## that displacement stays small against the period.
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
        fp = pss.factored_period()
        if fp.is_glm:
            raise NotImplementedError(
                "PAC.sampled_noise: the period map is '%s' (method %r); the "
                'seeded reverse pass exists for the linear multistep (gear, '
                'euler, trap) and stage (trbdf2 and other DIRKs, radau) maps, '
                "not for a multivalue GLM. Solve the PSS with method='radau' "
                "or 'gear'." % (fp.kind, getattr(pss.par, 'method', None)))
        stage = fp.is_stage
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

        ## components, rolled to the INJECTION index: step j's source enters
        ## at t_{j+1} (the reverse pass's own pairing) -- sampled one step
        ## early the held variance read 1.30 kT/C instead of 0.9987
        ## ⚠ WHERE THE SOURCE ENTERS: an LMM step's source enters at
        ## `t_{j+1}`, so `CY` is sampled at `x(t_{j+1})`; a stage method's
        ## enters at every stage abscissa `t_j + c_k h`, so at the STAGE
        ## states.  Measured on the sampler at 400 points: the end-of-step
        ## state for every stage read the held variance 0.867 kT/C (radau)
        ## against 1.000000 at the stage states, and it is the stage-state
        ## value that holds under refinement.
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
        ## the held variance is short by a factor that looked like a kernel
        ## constant (3.0 x the pure tail on a half-on switch) and was this.
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
        ## ⚠ ON A STAGED SOLVE THE SAMPLE'S ADJOINT IS BORDERED (2026-09-22,
        ## item 3 of the list after E7) -- the dual of the bordered forward
        ## solve, as `adjoint_sideband_row`'s: the operator is the TOTAL
        ## map's transpose, the sample is read at FIXED time (its costate
        ## carries `dtheta/dx_0^T Pk_fixed[k0]^T d`), and the source's own
        ## motion of the crossings enters as a third reverse pass carrying
        ## `-zeta_k W_k` at the event nodes, `zeta = Gt^-T (g_theta + a
        ## P_theta^T z)`.  Unbordered, the jitter sampler's held node had no
        ## path to the threshold's noise at all.
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
        `(N s,)` -- the injection times of a stage method's source."""
        tms = np.asarray(fp.times, dtype=float)
        out = []
        for j, st in enumerate(fp.steps):
            h = tms[j + 1] - tms[j]
            out.extend(tms[j] + st.c * h)
        return np.asarray(out, dtype=float)

    def _stage_states(self, pss, fp):
        """The stage states of one period, `N s` full-width vectors in the
        order of `_stage_times` -- one re-traversal of the converged orbit
        on a FRESH inner transient (the run's own is restored), cached per
        factored period."""
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
        """One reverse pass of a stage period map (`dirk` or `full`): returns
        the final costate and `(N s, m)` coupling vectors `h sum_i A_ik p_i`
        -- the sensitivity of the costate's functional to a unit source at
        stage `k` of step `j` is minus that (see `_sideband_forced`,
        whose loop this is).  `seed = (k0, v)` adds `v` to the costate after
        step `k0`'s update: the output at `t_{k0}` couples to the sources of
        earlier steps only."""
        m = pss.cir.n - 1
        N = len(fp.steps)
        lam = np.asarray(lam0, dtype=complex).copy()
        per = [None] * N
        for j in range(N - 1, -1, -1):
            st = fp.steps[j]
            lam, r = st.adjoint(lam)
            per[j] = [st.h * cp if cp is not None else np.zeros(m, dtype=complex)
                      for cp in (st.reach(r, k) for k in range(st.s))]
            if seed is not None:
                if isinstance(seed, dict):
                    if j in seed:
                        lam = lam + np.asarray(seed[j], dtype=complex)
                elif j == seed[0]:
                    lam = lam + seed[1]
        return lam, np.asarray([v for row in per for v in row], dtype=complex)

    @staticmethod
    def _psd_sqrt(Cs):
        """Symmetric PSD square roots of a stack `(..., n, n)`."""
        Cs = np.asarray(Cs, dtype=complex)
        Cs = 0.5 * (Cs + np.conj(np.swapaxes(Cs, -1, -2)))
        lam, U = np.linalg.eigh(Cs)
        return np.einsum('...ik,...k,...jk->...ij', U,
                         np.sqrt(np.clip(np.real(lam), 0.0, None)), U.conj())

    ## ⚠⚠⚠ THE NOTE THAT WAS HERE ACCUSED THIS ROUTE AND WAS WRONG.  A
    ## MONTE CARLO SETTLED IT THE OTHER WAY: this route is CORRECT and
    ## `orbital_correlation` is the one that fails on an asymmetric orbit.
    ## Direct SDE simulation of the variational system (trapezoidal, the
    ## calibrated `Var(i) = CY/(2h)` injection, phase projected out every
    ## step), sharing no Lyapunov solve and no modal sum:
    ##
    ##     a      |R| MONTE CARLO  |R| modal    |Lyap| proj   MC/Lyap
    ##     0.00   4.3629e-06       4.4515e-06   4.4437e-06    0.9818
    ##     0.30   2.9992e-04       3.6913e-06   2.9986e-04    1.0002
    ##
    ## `a = 0` is the CONTROL -- both routes agree there, so the MC had a
    ## known answer to hit, and it did (2 %).  At `a = 0.30` it lands on this
    ## route to 0.02 % and is 81x from the modal one.
    ##
    ## ⚠⚠ THE ARGUMENT THAT MISLED ME, RECORDED BECAUSE IT WAS PLAUSIBLE:
    ## `|lam2|` FALLS with asymmetry, so relaxation gets FASTER, so the
    ## transverse variance "should" shrink -- and the modal route did shrink
    ## while this one grew 67x.  That reasoning is WRONG: asymmetry changes
    ## the MODE SHAPES, so the noise projected onto the orbital direction
    ## grows, and the variance rises DESPITE the faster relaxation.  A
    ## physical argument is not a measurement.
    ##
    ## The original (refuted) note follows for the record:
    ## ⚠ SUPERSEDED: THIS ROUTE DEPARTS FROM PHYSICS
    ## ON AN ASYMMETRIC ORBIT, AND `orbital_correlation` DOES NOT.  van der
    ## Pol + `a u^2`, sweeping `a` (orbit asymmetry 0 -> 0.41):
    ##
    ##     a      |R| modal    |Lyap| proj  |K_orb raw|  |lam2|    amp
    ##     0.00   4.4515e-06   4.4437e-06   5.0034e-06   0.882521  2.000
    ##     0.05   4.4702e-06   4.7152e-06   4.7471e-06   0.881719  2.004
    ##     0.15   4.9505e-06   1.7371e-05   2.7846e-05   0.874886  2.034
    ##     0.30   3.6913e-06   2.9986e-04   9.5944e-04   0.844322  2.164
    ##
    ## `|lam2|` FALLS (0.883 -> 0.844), so amplitude relaxation gets FASTER
    ## and the transverse variance should get slightly SMALLER.  The modal
    ## sum does exactly that (4.45e-06 -> 3.69e-06).  This route grows 67x
    ## projected and 192x raw.  ⚠ The RAW covariance grows MORE than the
    ## projected one, so it is not the oblique projection -- it is this
    ## covariance.
    ##
    ## ⚠ FOUR EXPLANATIONS EXCLUDED BY MEASUREMENT, not by argument:
    ##   * harmonic truncation -- the disagreement is FLAT at 9.877e-01 from
    ##     `H = 4` to `H = 128`;
    ##   * the projection's tangent proxy -- `|cos(u, tangent)| = 1.000000`
    ##     at every asymmetry, against the Floquet phase mode;
    ##   * the modal decomposition -- `|lam1| = 1.000000`, `lam2` real and
    ##     well separated, one orbital mode, `p(T)-p(0) ~ 1e-14`,
    ##     `q^T C p = 1.000000`;
    ##   * a defect in `orbital_correlation` -- its two internal routes agree
    ##     to 3.4e-04 independently of `a`.
    ##
    ## ⚠⚠ WHAT THIS DOES **NOT** INVALIDATE.  Every use of this function as a
    ## reference in this file was on a SYMMETRIC orbit, where the two routes
    ## agree to 0.3 % -- including A9 step 3's three-way gate and the C^2
    ## biorthonormalisation defect it caught on 2026-09-07.  Those stand.
    ## ⚠ WHAT IS OPEN: which route is right is NOT settled.  The physical
    ## argument favours the modal one, but that is an argument, and this file
    ## does not close items on arguments.  A transient MONTE CARLO of the
    ## orbital fluctuation is the decisive third route and has not been run.
    ## Until then, treat this on a strongly asymmetric orbit as unvalidated.
    def oscillator_covariance(self, pss, samples=False):
        """
The state covariance of a FREE-RUNNING oscillator, split in two.

        Returns `(K_orb, d, info)`.  `K_orb` is the BOUNDED periodic
        (orbital) part of the covariance at `t = 0`; `d` is the growth per
        period along the orbit tangent, so

        ⚠ "BOUNDED" IS NOT "TRANSVERSE".  `K_orb` has the SECULAR growth
        removed and still contains the phase direction's bounded
        within-period variance.  Demir's orbital deviation `y` is the
        OBLIQUE projection `v_1^T y = 0`, so the transverse covariance is
        `Pi K_orb Pi^T` with `Pi = I - u v^T/(v^T u)` -- which is what
        `orbital_correlation`'s eq (23) sum equals (to 1e-4), and what
        `K_orb` itself does NOT equal (2-6 %, falling as 1/Q).  Read
        `K_orb` as the bounded part; project it if you want `R_yy(0)`.
        See `orbital_correlation`.

            K(t_0 + n T) = K_orb + n d u u^T

        exactly, for every integer `n`, with `u` the pair-space tangent
        scaled so its first block is `xdot(0)`.

        ⚠ WITH `samples=True` THE SPLIT MOVES WITH THE ORBIT, AND IT IS
        WORTH SAYING PRECISELY BECAUSE THE OBVIOUS READING IS WRONG.
        `info['orbital_samples'][j]` is `P(t_j)`, the solution started from
        `K_orb` at `t = 0`, and it satisfies

            K(t_j + n T) = P(t_j) + n d u_j u_j^T,   u_j = Phi(t_j, 0) u

        so `P` is periodic UP TO the growth -- `P(T) = P(0) + d u u^T`, not
        `P(T) = P(0)`.  The walk is along the orbit and the orbit turns, so
        the growth DIRECTION is the propagated tangent rather than a fixed
        `u`.  MEASURED: `P(T) - P(0)` matches `d u u^T` to 3.2e-15, and the
        full prediction holds to 2.5e-09 against a brute-force recursion
        run forty periods (9,600 steps) from `K = 0`.

        ⚠ THIS IS THE OBJECT `covariance` REFUSES TO RETURN, and the
        refusal was right: there is no periodic solution, so anything that
        returned one number would be hiding the physics.  `lambda_1 = 1`
        gives `lambda_1^2 = 1`, so `I - M kron M` is exactly singular --
        MEASURED here at `sigma_min` 2.3e-11 with the next singular value
        at 0.997, i.e. a null space that is cleanly ONE-DIMENSIONAL and
        spanned by `u kron u`, with left null `v kron v`.  So it borders
        exactly as the PPV and the deflated PAC solve do, and the border is
        the pair the rest of this class already computes.

        ⚠ THE SPLIT IS NOT A NUMERICAL DEVICE, IT IS THE ANSWER.  Demir
        2002: an oscillator's noise is STATIONARY, not cyclostationary,
        because "noisy autonomous systems cannot provide a perfect time
        reference".  `K_orb` is the part a designer can read as an
        amplitude/orbital noise -- it settles, it is periodic, it is
        finite.  `n d u u^T` is the random walk ALONG the orbit, which
        never settles and which no periodic object can hold.  Reporting
        only their sum at some finite time is what makes an oscillator
        covariance look divergent and useless; reporting the parts makes
        both usable.

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
        `ppv()` normalises on the FIRST BLOCK, `v[:m] . xdot = 1`, which is
        the normalisation a state perturbation entering the first block
        sees -- an injected current, and what every other shipped path
        does.  The FULL PAIR contraction is a different number: 0.663 on
        van der Pol, so `(v . u)^-2 = 2.28`.  That exact mistake produced a
        2.31x discrepancy that was chased as a code defect for a while; it
        is why `d` is written with the pair inner product spelled out.

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
        Neither touches the other's machinery.

        ⚠ AND THE TWO ANCHORS BEHIND THEM ARE ALSO INDEPENDENT, which is
        the property that was missing when a 2x error survived a 0.9965
        agreement.  `covariance`'s injection is anchored to `kT/C`;
        `diffusion_constant` is anchored to a nonlinear Monte Carlo reading
        phase from zero crossings.  `info['c_from_growth']` against
        `diffusion_constant` therefore closes a loop between two separately
        anchored quantities rather than reproducing one of them.

        ⚠ COST: the bordered solve has `(2m)^2 + 1` unknowns and is dense,
        so it is `O(m^4)` like `covariance`.  Small circuits only.  The
        closed form for `d` is cheap; pass `samples=False` and read
        `info['c_from_growth']` if the orbital part is not wanted.
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
        ## noise-driven motion in the injection (2026-09-22) -- the same
        ## `_event_closure` as `covariance`, whose `u`, `v` below are the
        ## total map's already
        _bordered = self._event_closure(pss, As, Qs, M, m, n)
        if _bordered is not None:
            M, K1, _samples_unused, _pieces_unused = _bordered

        v, pinfo = pss.ppv()
        v = np.asarray(v, dtype=float).ravel()
        u = np.asarray(pinfo['tangent_pair'], dtype=float).ravel()
        xdot = np.asarray(pinfo['xdot'], dtype=float).ravel()
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
            for A, Q in zip(As, Qs):
                K = A @ K @ A.T + Q
                uj = A @ uj
                orb.append(0.5 * (K + K.T))
                grw.append(d * np.outer(uj, uj))
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
        METHOD: measured 0.1611 / 0.0173 / 0.0017 / 0.0002 as the TANK's noise
        falls 1e-6 -> 1e-9 against a fixed 1e-6 at the buffers, because the
        phase direction is exactly what tank noise drives.  So a fixture whose
        oscillator is quiet will show the projection doing nothing and teach
        you the wrong lesson; one with a noisy tank shows 16 %.  It also falls
        as 1/Q, so a high-Q fixture hides it too.
        ⚠ `P - (t/T) G` is the SUPERSEDED prescription and is also wrong; see
        `orbital_correlation`, which gates the projection three ways.

        ⚠ `u_j` COMES BACK FROM `growth_samples[j]`, WHICH IS A RANK-ONE
        MATRIX `d u_j u_j^T`, not a vector -- its leading eigenpair gives
        `sqrt(d) u_j`, and `Pi` is invariant to that scale (and to `v`'s), so
        taking `v` from `ppv()` in a separate call is safe here.  Everything
        is sliced `[:m, :m]` out of PAIR space.

        ⚠ THE SLOPE IS A LOCAL QUADRATIC FIT, NOT A TWO-POINT DIFFERENCE, and
        that is a correction rather than a preference.  A straddling
        difference across a threshold crossing samples whichever pair of grid
        points brackets it, and the crossing sits at a different fraction of a
        step on every grid: measured -0.720 / -0.077 / +2.028 % under
        refinement, changing SIGN, while the period and swing of the same runs
        converged cleanly at first order.  The quadratic fit reads
        1.527458 / 1.527793 / 1.528023 at 240 / 480 / 960 points -- flat.  The
        straddling value on a 240-point grid was 2.8 % low, and since every
        quantity here goes as `1/s^2` that inflated the variance by 5.6 % and
        was briefly blamed on the other side's integrator.

        MEASURED 2026-09-17 on a van der Pol tank driving three tanh buffers
        (`tau_buf = RC = 0.5`), noise at every node: `2 sigma_t^2` reads
        1.055708e-06 / 1.055112e-06 / 1.062936e-06 at 240 / 480 / 960 points,
        i.e. FLAT to a few parts per thousand.  ⚠ An earlier record of this
        said it CONVERGED (1.023870e-06 -> 1.055112e-06 -> 1.062936e-06, "the
        increments shrinking ~4x").  That was the slope error above shrinking
        with the grid, not the covariance converging -- with the slope taken
        at the requested instant the grid dependence is essentially gone.  A
        Monte Carlo over EIGHTY
        seed-runs -- a noisy transient with no PSS, no adjoint and no Lyapunov
        solve in it -- gives

            MC / analysis = 1.0066 +/- 0.0102   (0.64 sigma from 1.000)

        over 124 seed-runs across three grids, ONCE the comparison is made
        between the same orbit on both sides.

        ⚠⚠ AN EARLIER RECORD OF THIS CLAIMED A 4-SIGMA DEFECT IN THIS METHOD
        AND IT WAS WRONG -- THE FAULT WAS IN THE COMPARISON, NOT HERE.  An
        80-seed campaign reported `1.0571 +/- 0.0141`, "under-predicts by
        5.7 %, grid-independent, therefore on the analysis side".  It was
        neither.  `sigma_t = sqrt(var)/slew`, so the comparison goes as
        `1/slew^2` -- and the Monte Carlo finds its crossings on an EULER
        orbit while this method divides by the slope of the PSS's GEAR orbit.
        At finite `h` those orbits differ: Euler's period runs 0.39 % short
        and its slope at the crossing is low by 2.717 / 1.361 / 0.686 % at
        240 / 480 / 960 points, HALVING as first-order convergence requires.
        Squared, that predicts ratios of 1.0566 / 1.0278 / 1.0139.

        MEASURED AGAINST THAT PREDICTION (the 960 value pinned before its
        seeds ran): raw 1.0493 / 1.0513 / 1.0165, and dividing by each grid's
        own independently measured slew correction collapses them onto
        0.9931 +/- 0.0153, 1.0229 +/- 0.0160, 1.0025 +/- 0.0261 -- pooled
        1.0066 +/- 0.0102, with a residual trend of +0.0047 per doubling
        against a per-grid scatter of 0.019.

        ⚠ THE ERROR THAT PRODUCED THE FALSE ALARM IS WORTH MORE THAN THE
        NUMBER: "grid-independent" was asserted from TWO points 0.99 sigma
        apart.  Two noisy points cannot tell FLAT from HALVING, and halving is
        what it was doing.  An absence of evidence was read as evidence of a
        property, and a structural conclusion ("therefore the analysis's") was
        built on it.

        A factor of 2 in the variance remains excluded at more than 15 sigma.
        ⚠ When validating this against a transient, run the Monte Carlo on the
        SAME integrator as the PSS, or divide by the slope of the orbit the
        Monte Carlo actually runs on -- otherwise the mismatch enters squared.

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
        """
        import warnings as _warnings
        self._check_circuit(pss)
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
        ## instead reproduces the straddling value (1.485472 against a
        ## converged 1.5275 on a 240-point grid, 2.8 % low) and every quantity
        ## here goes as 1/s^2.
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
        ## ⚠ `samples[j]` IS node j (2026-09-20, measured).  This list used to
        ## be `[v0] + samples`, which paired node j's covariance with the
        ## phase vector of node j - 1 -- a one-node shift, first order in the
        ## step: the obliquely-projected transverse cycle mean read 8.0e-2 /
        ## 4.0e-2 / 1.7e-2 against its own N = 3200 value at N = 200 / 400 /
        ## 800 (halving), and 3.3e-4 / 5.1e-4 / 1.8e-4 with the list
        ## unshifted.  The same prepend sat in three tests, and it was read
        ## as "the Lyapunov route is first order" (5129485) -- it was this.
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
        returns the NON-NULL directions, so `sum cw[k,k'] u_k u_k'^H` reproduces
        only the part of `K_orb` that lives on them.  How much that is depends
        entirely on WHERE THE NOISE ENTERS -- measured on `_osc_with_ladder`'s
        circuit at `nslow = 4`, moving one current source and changing nothing
        else::

            injected at            ||K_orb||    reconstruction residual
            the oscillator node    2.70e-05     1.80e-03   (0.18%)
            a SLOW ladder node     3.94e-01     3.56e-01   (36%)
            a FAST ladder node     6.87e+02     9.996e-01  (99.96%)
            a faster one           3.33e+03     9.999e-01  (99.99%)

        **When the injection lands in a fast branch the non-null modes capture
        essentially NOTHING of the covariance.**  The annihilated modes are
        killed by the period map, so they enter the stationary covariance only
        through the `j = 0` term -- but that term is not small when the noise
        is injected there, and THAT IS WHERE DEVICE NOISE ACTUALLY IS: every
        resistor in a bias or tuning network.

        ⚠ SO A MODAL ORBITAL SPECTRUM BUILT ON THIS BASIS IS COMPLETE ONLY FOR
        NOISE THAT ENTERS THE SLOW SUBSPACE, and the suite's own gate on this
        (`rel < 1e-2`) holds because its fixture injects at the oscillator
        node.  That is a property of the fixture, not of the method.

        ⚠ AND THE RESIDUAL IS A DETECTOR, NOT A TRUNCATION BOUND.  It catches a
        DROPPED NON-NULL MODE well -- which is what the note below claims for
        it -- but it SATURATES at the floor above, so it cannot certify a
        truncation below whatever the null modes carry, however many modes are
        kept.  Independently reproduced by a peer session on a different
        oscillator with a different `K_orb` route (69% there, mechanism
        identical, magnitude not transferable).

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
        eq (22), and the OUTPUT spectrum then needs a layer that can carry
        contributions asymmetric about the carrier. Those are not built.
        """
        K_orb, _d, _info = self.oscillator_covariance(pss)
        K = np.asarray(K_orb, dtype=float)
        n = K.shape[0]
        modes = pss.floquet_modes(pss, nmodes=(n if nmodes is None
                                               else int(nmodes)))
        V = np.column_stack([m['v0'] for m in modes])
        cw = V.conj().T @ K @ V
        return cw, modes, K

    ORBITAL_HARMONICS = 32

    ## Half-wave asymmetry above which `orbital_correlation` is known to be
    ## wrong.  MEASURED (below); 0.02 is a decade inside the smallest
    ## asymmetry at which the error was already visible.
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

        ⚠⚠ THIS WARNING WAS WRITTEN FOR A DEFECT THAT IS NOW FIXED, and kept
        for the residual.  On 2026-09-07 `orbital_correlation` read 81x LOW
        against a Monte Carlo on van der Pol + `a u^2` at asymmetry 0.41.
        The cause was in `floquet_modes`: the replayed adjoint is `C^T q`,
        not `q`, and was used untransformed -- invisible on a unit-reactance
        symmetric orbit, catastrophic off-axis where the two adjoints are
        nearly parallel.  With the `C^-T` transform applied, against the
        same Monte-Carlo-validated Lyapunov reference at `a = 0.30`:

            npts    eq22 / Lyapunov
             400       1.0595
             800       1.0300
            1600       1.0151

        halving per doubling -- an `O(h)` DISCRETISATION residual of the
        adjoint replay, converging to 1, not a defect.  At `a = 0` it is
        1.0001.  So this warns that the residual is grid-limited on such an
        orbit and says how to shrink it; it no longer says the answer is
        wrong, because it is not.
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
        integrates `B B^T` as a two-sided intensity.  Consistent with the
        `kT/C`-calibrated Monte Carlo injection `Var(i) = CY/(2h)` in this
        file's record, and confirmed here by three routes agreeing.

        ⚠⚠ GATED THREE WAYS, because a modal sum transcribed from an image
        of an equation is exactly the object this record distrusts.  On
        van der Pol under gear: (i) this sum against `R_yy(0)` evaluated
        from its DEFINITION as a 1-D Lyapunov integral along the orbital
        mode, no Fourier machinery — agree to 3.5e-4; (ii) both against
        the CYCLE-MEAN transverse part of `oscillator_covariance`'s
        samples (`P(t_j) - (t_j/T)·growth_samples[j]`), which shares no
        machinery with either — magnitude to < 1e-3.  That third route is
        what found the state-block scale defect in `floquet_modes`.

        ⚠ THE REFERENCE IS THE CYCLE MEAN, NOT `K_orb(0)`.  Lemma 3.5's
        `R∞_yy` depends on `τ` only — the stationary part.  At `t = 0`
        van der Pol's amplitude direction is pure-v while this is
        isotropic, which is a rotating radial direction averaged over a
        cycle, not a disagreement.

        ✅ THE 2-3 % SHAPE RESIDUAL WAS THE REFERENCE, NOT THIS SUM -- closed
        2026-09-04.  Subtracting only the SECULAR growth `(t/T) d u u^T`
        from the Lyapunov samples leaves the phase direction's BOUNDED
        within-period variance, which eq (22)'s `l >= 2` sum correctly
        excludes.  Demir's `y` is defined by the OBLIQUE projection
        `v_1^T y = 0`; project the samples with `Pi = I - u v^T/(v^T u)`
        and the three-way agreement is 5.9e-4 / 6.0e-4 / 3.1e-4 / 1.5e-4
        at Q = 4 / 8 / 16 / 32 (euler-plain, n = m), improving with
        refinement -- quadrature.  The old residual fell as 1/Q_lambda
        (5.4 / 2.4 / 1.2 / 0.6 %), which is orbital variance ~ Q against a
        constant phase-bounded part: the same fact, seen from the sweep.
        ⚠ A candidate recorded earlier -- the phase-orbital CORRELATION's
        tau = 0 value -- was DISPROVED from eq (18a) before it was built:
        at tau = 0 its brace is {1 - 1} = 0 identically, and (23) states
        R_yy(0) = sum C_lhj alone.  Named so nobody rebuilds it.
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
        ## blew up as 1/|mu|^2 (146x .. 2449x radau's R on a 3:1 gear grid)
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

        ⚠⚠ MEASURED 2026-09-14 -- "negligible" is a property of THOSE
        circuits, and the over-statement can be large.  Against pnoise (the
        total linear sideband noise, on one absolute scale with the certified
        Lorentzian), `(up+lo)/(4 (S_ph + S_orb))` at 10 f_amp on van der Pol
        with C = 4, Q = 8 and an `a u^2` asymmetry:

            half-wave asymmetry   0      0.017   0.033   0.067   0.100
            R                     0.999  0.997   0.972   0.691   0.309

        So on a SYMMETRIC orbit the amplitude is right to 0.1 % (also at
        C = 1 and Q = 50 -- the first external check this spectrum had), and
        on an asymmetric one the sum over-states by up to 3.2x (5 dB) with no
        grid dependence.  A Monte Carlo of the SDE (64 oscillators x 4000
        periods, trapezoidal, `Var(i) = PSD/(2h)`) settles which side is
        right: at a = 0.30, MC/pnoise = 1.011 and MC/(S_ph + S_orb) = 0.313,
        with the a = 0 control reading 1.009 / 1.008.  pnoise is the total.

        ⚠⚠ THE DROPPED CROSS TERM IS THE CAUSE -- BUT ONLY WITH EVERY
        HARMONIC KEPT (2026-09-15; the reverse reading of 2026-09-14 held only
        for the truncated form).  `S_corr` built from eq (92) on this fixture
        is ~1e-8 of the total: that form keeps only the PPV's DC harmonic AT
        THE NOISE SOURCE's row, and an ideal tank inductor shorts that node at
        DC (`vbar` = [-2.3e-6, -0.114]).  The full-harmonic correlation is
        -1.1 to -2.4x this spectrum at a = 0.30, and `modal_spectrum`'s three
        terms sum to pnoise (1.006 at 10 f_amp, grid error that halves with
        the grid).  The orbital mode's AM share at the output,
        `sin^2 arg(U_{l,1}/U_{0,1})` = 0.307, is the flat factor below.
        The true total is BELOW even the phase Lorentzian alone (pnoise/S_ph
        = 0.61 at a = 0.30): the over-statement is in the decomposition's
        frequency-independent terms above f_amp, not in a missing
        correction.  LOCALISED: the PHASE half is the Lorentzian's
        frequency-independent PPV -- with `c(f)` from `frequency_aware_ppv`
        pnoise's PM content matches it to <= 2.3 % at 0.3-10 f_amp -- and
        the ORBITAL half over-states by a factor FLAT in offset (AM content
        0.317 of this spectrum at a = 0.30, at every offset), which is open.
        Traversa & Bonani's own Figs 1-2 show the same limit on their
        amplitude-phase-coupled test oscillator (theory above the exact
        spectrum at high frequency, growing with the coupling).  A warning
        fires above `ORBITAL_ASYMMETRY_LIMIT`.

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
        """
        try:
            _asym = self._orbit_asymmetry(pss)
        except Exception:
            _asym = 0.0
        if _asym > self.ORBITAL_ASYMMETRY_LIMIT:
            ## ⚠ NOT the grid residual `_warn_if_orbit_is_asymmetric` names:
            ## measured against pnoise (2026-09-14), the sum this spectrum is
            ## meant for over-states the TOTAL on an asymmetric orbit, and no
            ## refinement changes it -- see the docstring.
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
        ## tails of the others (2026-09-14).  The line weight at `j f0` is
        ## `W_j = sum_{l,h} Re(row C_lhj row)`; where it is zero (a symmetric
        ## orbit's even harmonics: the modes' own Fourier content vanishes)
        ## what this would return is the neighbouring lines' Lorentzian tails,
        ## which a Monte Carlo put 3.2x LOW at 2 f0 on van der Pol C=4 Q=8
        ## (pnoise agreed with it to 1 %).  ⚠ NOT caught: DC, where a small
        ## line can exist and the model read ~100x HIGH (the tank inductor
        ## shorts the node, which Lorentzian tails do not know), and 2 f0 on
        ## an asymmetric orbit (0.40) -- away from the fundamental use
        ## `pnoise`.
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
        transfer, which sum to the total.  E6, built 2026-09-15.

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
        3.2x (Monte-Carlo-confirmed).  The missing piece IS the correlation --
        but not in the form Traversa & Bonani keep: their eq (92) retains only
        its DC harmonic, ~1e-8 of the total on van der Pol (the tank inductor
        shorts the source node at DC).  With every harmonic kept it is -1.1 to
        -2.4x the orbital term there.  It removes the DC-PPV phase excess above
        f_amp AND the orbital mode's PM projection at the output -- the orbital
        line's AM share `sin^2 arg(U_{l,1}/U_{0,1})` (0.307 at a = 0.30) is the
        flat factor `orbital_spectrum` was measured to over-state by.
        Measured on van der Pol C=4, Q=8, 400 points per period (H = 8, 16
        sidebands), `total / (pnoise/2)`:

            a     harmonic   +1      +3      +10     -3      -10  f_amp
            0.00  1          1.0006  1.0006  1.0006  1.0006  1.0006
            0.00  2                  1.0010  1.0010  1.0010  1.0009
            0.30  1          1.054   1.014   1.006   1.008   1.006
            0.30  2                  1.014   1.005   1.008   1.005

        ⚠ THE a = 0.30 EXCESS IS GRID ERROR, NOT THE MODEL: on 800 points it
        reads 1.026 / 1.006 / 1.002 at +1/+3/+10 f_amp (it halves -- the O(h)
        of the modes), and the a = 0 control reads 1.00013; H 8 -> 12 and
        sidebands 16 -> 24 move the fourth digit.  ⚠ Because the correlation
        cancels most of the other two, a few percent of error in any part is
        AMPLIFIED in the total -- which is why the three are computed together
        here rather than the correlation being offered as an add-on to
        `oscillator_spectrum + orbital_spectrum`: those line-shape spectra keep
        only the resonant term of each line (2.5 % short at 10 f_amp even on a
        symmetric orbit), which is harmless alone and not under cancellation.
        The modal sum also reproduces pnoise's upper/lower sideband asymmetry,
        which the two-term sum cannot.

        Near the carrier `phase` IS the library Lorentzian (1.00022 of
        `oscillator_spectrum(frequency_aware=False)` from 0 to 3 linewidths,
        symmetric orbit) and `correlation` is ~1e-6 of it.  Above f_amp
        `total` agrees with `pnoise`; within the linewidth pnoise has no
        meaning and this is the route.

        ⚠ Stationary WHITE sources, free-running oscillators, and the dense
        `floquet_modes` only (inherited).  `harmonic >= 1`: harmonic 0 was
        never measured.  `H` defaults to `ORBITAL_HARMONICS` (capped by the
        grid), `sidebands` to `2 H`.  `output` follows `orbital_spectrum`.
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
        ## `_phase_mode_split` (the 1e-6 window refused every gear solve on a
        ## non-uniform grid)
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
            ## ⚠ `_period_dft`, not an index DFT: measured 8-13 % off and NOT
            ## converging on a 3:1 grid before (see `PSS._period_quadrature`)
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
        a window on `|lam| - 1` (2026-09-20, Andreas: gear as a first-class
        choice on non-uniform grids).

        On a uniform grid the phase multiplier is 1 to rounding.  On a grid
        whose step varies, a multistep or trapezoidal solve loses time-
        translation symmetry and the multiplier leaves the circle at O(h^2)
        -- measured 1 - 5.1e-05 (gear) and 1 + 4.7e-05 (trap) at 400 points
        on a 3:1 grid; radau keeps it at 1 + 1e-11 -- and a 1e-6 window
        refused every gear solve there.  The right eigenvector of the phase
        mode is the tangent `xdot(0)` (`C xdot = -i(x)` for the autonomous
        circuit), which no orbital mode shares, so alignment picks it on any
        grid; its exponent is then forced to 0 exactly, as the consumers
        already do.  The departure is WARNED with its size when it exceeds
        rounding, and the split is REFUSED when a second multiplier lies
        within ten times that departure of the circle with any alignment --
        the case a window ever protected against.
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

        ⚠ ITS SCALE WAS WRONG BY 2x AND IS NOW FIXED -- kept because the
        way it survived is the instructive part.  `diffusion_constant` used
        the full `CY` while `covariance` used `CY/2`: two functions in this
        class disagreeing about whether `CY` is one- or two-sided.  It was
        validated against a Monte Carlo injecting `Var(i) = CY/h` per step
        and agreed to 0.9965 -- because that Monte Carlo carried the SAME
        hot convention.  A measurement built on the assumption under test
        cannot test it.

        SETTLED AGAINST `kT/C`, which is external to both: an injection of
        `Var(i) = CY/h` reproduces 1.92x `kT/C` over ten independent runs
        (1.75-2.04).  With `CY/2` throughout, `diffusion_constant` gives
        7.9516e-08 against a correctly scaled Monte Carlo at 7.7083e-08 --
        ratio 1.0316, inside that measurement's 4.1% uncertainty.

        ⚠ AND A SECOND DISCREPANCY WAS NOT A CODE DEFECT AT ALL.  Two Monte
        Carlo routes disagreed by 2.31x, which looked like a third error.
        It was in the DIAGNOSTIC: `ppv()` normalises on the FIRST BLOCK
        (`v[:m] . xdot = 1`), which is right for a perturbation entering
        the first block -- an injected current, and what every shipped path
        does -- but wrong for contracting against a full PAIR deviation,
        where the factor is `1/(v . u_pair) = 1.508`.  Correcting it turned
        a 2.13 variance ratio into 1.07.  The sign difference alongside it
        is a convention, not an error: a later zero crossing means DELAYED,
        while projecting onto the tangent makes positive mean ADVANCED.
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

        ⚠ EVERY PERIOD INTEGRAL HERE USED `h = diff(times)` -- the LEFT
        RECTANGLE rule (2026-09-20).  On a uniform periodic grid that is the
        trapezoid rule, spectrally accurate; on a NON-UNIFORM grid it is
        FIRST order by itself (its error is (1/2) integral h'(t) y(t) dt, not
        zero).  Measured on a solved-history (gear) run of an index-2
        oscillator on a smooth grid, `c` against radau: -1.9e-3 / -9.3e-4 /
        -4.6e-4 at N = 200 / 400 / 800 with the rectangle weights and
        +4.0e-5 / +3.6e-5 / +1.3e-5 with these, from the same samples.  On a
        uniform grid `0.5 h + 0.5 h == h` exactly, so every uniform-grid
        number is bit-identical to before.  The one-step kinds' replays are
        uniform-grid replays (`factored_period_stage`), so only the
        solved-history kind ever paid this.  ⚠ AND THE TRAPEZOID IS ITSELF A
        SECOND-ORDER CAP on a smoothly varying grid (2026-09-21): with `pss`
        given, a non-uniform grid and no landed events, these are the
        periodic cubic-spline weights of `periodic_spline_weights` (radau's
        `c` on a 1 + 0.5 sin grid 1.3e-5 -> 2.4e-10 at N = 200); under
        landed events the spline breaks at their nodes (2026-09-21, it used
        to keep the trapezoid), and a uniform grid is unchanged."""
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
        value `diffusion_constant` used to return silently.  `phase_psd`
        reads it at the carrier for the Lorentzian CORNER, which is a
        white-noise construct whatever the source's colour; the spectrum
        itself comes from `coloured_diffusion_resolved`.
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
        ## covariance and `samples` is `C^T v_1`; contracting it here made
        ## `c` wrong by `C^2` on the differential rows and exactly zero on
        ## the algebraic ones.  See `_equation_row_ppv`.
        S = np.asarray(info['samples_eq'])[:, :m]
        tms = np.asarray(info['times'], dtype=float)
        ## the samples' own orbit, not `pss.period` -- see `ppv()`'s 'period'
        T = float(info['period'])
        h = self._period_weights(tms, S.shape[0], T, pss)
        cy = self._cy_reduced(pss, float(w))
        ## ⚠ A NOISE SOURCE ON AN INDEX-2 CONSTRAINT GIVES c = 0, SILENTLY
        ## (2026-09-21): a voltage noise in series with a DC source inside a
        ## capacitor loop perturbs an algebraic constraint -- a DIFFERENTIATED
        ## input, whose response is a charge jump the PPV projection cannot
        ## represent -- and `c` came back exactly 0 for every method on an
        ## index-2 van der Pol (the same circuit's current noise at the node
        ## gives 3.2e-9).  Named here once: the PPV's algebraic fallback
        ## fired (index >= 2) and `CY` has power on an algebraic row.
        ## (the gate is the index-2 condition itself -- `G[A,Z]` singular at
        ## the orbit point, the test the PPV's algebraic fallback makes --
        ## computed here so every kind is covered, not only solved-history)
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
        ## USES.  `CY` is a one-sided density (a resistor's `4kT/R`), and
        ## these two functions disagreed about it until a Monte Carlo was
        ## run against `kT/C`: an injection of `Var(i) = CY/h` per step
        ## reproduces `1.92x kT/C` over ten independent runs (1.75-2.04),
        ## so that convention carries TWICE the physical noise power.
        ## `covariance` was already right; this was not, and its agreement
        ## with a Monte Carlo built on the SAME hot convention is exactly
        ## why the error survived.
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
        answer is not, so nothing downstream would look wrong.  The
        measured separation is 22 ORDERS on van der Pol -- `c = 7.95e-08`
        against `Gamma = 1.9e-29` -- so the two functionals are not close
        approximations of each other and cannot be substituted.

        ⚠ TWO INDEPENDENT MECHANISMS FORCE `vbar` TO ZERO, AND ONLY ONE OF
        THEM IS THE ONE DESIGNERS KNOW.  Measured on an LC oscillator,
        sweeping an even term `a (u^2 - 2)` in the nonlinearity and a
        series tank resistance `Rs`:

            a      Rs      Gamma/c
            0.00   0.00    2.4e-22
            0.00   0.20    9.7e-23
            0.25   0.00    4.9e-23
            0.25   0.05    2.1e-04
            0.25   0.20    4.1e-03

        NEITHER ASYMMETRY ALONE NOR LOSS ALONE UPCONVERTS.  `c` is
        7.9e-08 to 1.2e-07 in every row, so the quadratic functional
        cannot produce that pattern.

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
        the pair cannot drift the way `diffusion_constant` and
        `covariance` once did over exactly that factor.
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

        ⚠ THIS IS THE OBJECT `c + Gamma(f)` STOOD IN FOR, and the stand-in
        is wrong in two ways that the fixture could not show: `c` reads
        `CY` at ONE frequency (`2 pi / T`) as if it held at every harmonic,
        and `Gamma` is exactly the `l = 0` term of this sum, so `c + Gamma`
        counts `l = 0` twice.  Neither was visible on van der Pol, whose
        PPV at the tank node averages to zero (`|V_0|/|V_1| = 5e-13`: the
        inductor shorts the node at DC, so no core can bias it) -- the
        fixture shared the claim's assumption, failure shape 0b.

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
        a white source it is `c` exactly; the earlier `c + Gamma(f)` form
        counted the `l = 0` term twice and is retired.

        ⚠ THE CONVENTION IS PINNED BY `oscillator_spectrum`, NOT ARGUED.
        `lorentzian`'s far skirt is `i^2 f_0^2 c / f^2` exactly, and that
        object was gated by power conservation to 1.000000.  So this
        expression is the same quantity its tail already reports, with the
        coloured term added -- no second convention is introduced, which
        is the only reason a `S_phi` is shipped here at all after a
        one-sided/two-sided error cost this class a factor of two.

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
        ## ⚠ THE CORNER IS THE WHITE LORENTZIAN'S, read at the carrier as it
        ## always was.  For a coloured source `f_h = pi i^2 f0^2 c` is not a
        ## lineshape parameter at all -- there is no Lorentzian -- and
        ## taking the folded value nearest the carrier instead put a 1/f
        ## source's corner ABOVE the offsets, in front of the power bound
        ## below, which is the floor that actually binds for colour.
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
        ## statement for a 1/f input and reduce it to
        ## `df_c >= eps f0 sqrt(2 f_1f)`; the form here needs no
        ## assumption about the source's colour, and REPRODUCES their
        ## worked example exactly -- 100.000 Hz against their ">= 100 Hz"
        ## at `eps^2 = 1e-19`, `f0 = 1 GHz`, `f_1f = 50 kHz`.
        ##
        ## ⚠ THE LORENTZIAN CORNER ABOVE DOES NOT CATCH THIS.  It is built
        ## from `c` alone, so it knows nothing about a `Gamma(f)` that
        ## grows as the offset falls.  MEASURED on this class's own
        ## flicker fixture: the power bound bites at 2.5e-06 Hz while the
        ## Lorentzian corner sits at 8.2e-09 Hz -- 306x too permissive,
        ## and the swept spectrum was carrying 3.10x unit power at the
        ## bottom of the range before this check existed.
        ##
        ## ⚠ IT IS A LOWER BOUND ON THE BREAKDOWN, NOT THE BREAKDOWN.
        ## Passing it is not a guarantee: on Vanassche's own example the
        ## observed flattening sits at ~300 Hz, 3x the bound.  So this
        ## refuses what is definitely invalid and admits a band that is
        ## already suspect -- deliberately, because refusing at 3x would
        ## be fitting a threshold to one example.
        ## ⚠ AND THE DERIVATION HAS A PRECONDITION THE BOUND DOES NOT
        ## STATE, so it is checked rather than assumed.  The box argument
        ## is `2 df S(df) <= integral_{-df}^{+df} S <= 1`, and the FIRST
        ## inequality needs `S(f) >= S(df)` for every `|f| <= df` -- the
        ## spectrum must not dip below its edge value anywhere further in.
        ## True of a monotone skirt; TRUE of the flattened near-carrier
        ## shape; true even with a spur, which ADDS power inside rather
        ## than creating a dip.
        ##
        ## ⚠ FALSE FOR A LOCKED PLL, whose phase-noise transfer function
        ## is HIGH-PASS: the spectrum is SUPPRESSED at DC and rises to the
        ## free-running level beyond the loop bandwidth, so it dips below
        ## its edge value everywhere inside.  The bound is not thereby
        ## shown to be violated there -- total power is still 1 -- it is
        ## NO LONGER DERIVED, and a floor that is not derived cannot be
        ## used as one.  Unreachable today because this method refuses a
        ## driven circuit, and squarely in the way of the driven-oscillator
        ## work, which is why it is a check and not a comment.
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
        two-regime approximations exist — which is one more reason this
        module supports white sources only.

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
        `1/r²`, so `S·r²` is flat and a band mean IS a point value — that is
        the case every gate in this tree was written on, and it makes the
        distinction invisible.  It is NOT general: MEASURED 2026-09-09, a
        source behind a slow RC node has an in-band spectrum that is not
        `1/r²` at all (its `k = 0` term is filtered at the RC corner while the
        `k >= 1` terms are not, and their mix moves across the band), and its
        slow/core ratio swings 1.16 -> 0.87 across `0.08 … 0.15 f0` — so a
        band mean and a point value differ by ~4 % there, which is larger
        than most of the agreements this file asserts.  A comparison that
        takes a band mean on one side and a point value on the other is then
        measuring the convention, not the physics.

        ⚠ So: call this before comparing a measured band-averaged number
        against a computed point value, or vice versa.  A spread near 1
        licenses the shortcut; anything else says put both sides on the same
        footing.  `quantity` selects the surface (`'pnoise'`, `'S_pm'`,
        `'S_am'`, `'oscillator_spectrum'`); `**kw` is forwarded to it.
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
        The fold is a FREQUENCY-CONVERSION computation, and for a driven circuit
        -- a mixer, the diode-mixer fold case -- that is complete.  For an
        AUTONOMOUS oscillator it is structurally incomplete, and what it omits is
        exactly the near-carrier phase-noise skirt this method returns.  Rizzoli,
        Mastri & Masotti (IEEE MTT 42-807, 1994) state it directly: frequency
        conversion alone is insufficient for autonomous circuits, because the
        noise-induced FREQUENCY MODULATION OF THE CARRIER at low offsets is not a
        frequency-conversion effect (verified at the source 2026-09-08: p. 807,
        Introduction, verbatim "frequency-conversion techniques alone are not
        sufficient to solve the noise analysis problem for general autonomous
        circuits (oscillators). An important further aspect that must be taken
        into account is the noise-induced frequency modulation of the carrier
        taking place at low frequency offsets, which is not a
        frequency-conversion effect").  Their Section III (p. 810) NAMES the
        two stacks: CONVERSION noise, power exchanged among the sidebands of
        the unperturbed steady state, "invariably raises as 1/f for f -> 0,
        which is not consistent with the measured behavior"; MODULATION noise,
        "a jitter of the oscillatory steady state", proportional to noise power
        over f^2 so the PSD "raises as 1/f^3 for f -> 0 in agreement with the
        measured performance"; and the two DECOUPLE exactly at the steady state
        (M_BH = M_HB = 0).  They also say the two are "usually nearly equal" in
        an INTERMEDIATE offset band, "so that (20) and (21) are
        interchangeable" -- a cross-stack agreement test this tree does not yet
        have (recorded in the roadmap, not built).  ⚠ Their construction is
        harmonic balance; what transfers is the classification, the two slopes
        and the interchangeability, none of which need HB.  Diagnostic value:
        a FLAT PSD near the carrier is neither slope -- it is the Phi(T) - I
        singularity, not the conversion model being the wrong physics.

        So the two stacks -- the Floquet/PPV one (`ppv`, `diffusion_constant`,
        this method) and the sideband fold (`pnoise`) -- ARE NOT TWO
        IMPLEMENTATIONS OF ONE QUANTITY, and unifying them is not a
        simplification waiting to be made.  ⚠ THE HAZARD IS THAT THE WRONG ONE
        STILL RETURNS A NUMBER: deriving oscillator phase noise from the fold
        alone would produce a spectrum -- the conversion terms are real and
        non-zero -- just one missing the dominant contribution near the carrier.
        A plausible wrong answer, which is the failure shape this whole area
        keeps generating.  That is the completeness argument for the split; the
        efficiency argument (Floquet is cheaper) is the weaker one and was for a
        long time the only one written down.

        Returns `(S_v, L_dBc)`.  ⚠ `S_v` is the Lorentzian lineshape scaled by
        `|X_1|^2 = A^2/4`, the carrier PHASOR's square -- which is HALF the
        carrier power `A^2/2` a one-sided PSD carries, so `S_v` is exactly
        0.5000x a one-sided PSD of the output voltage (measured against
        a reference simulator at every offset over four decades, 2026-09-05).  `L_dBc`
        is unaffected, `|X_1|^2` dividing out of the ratio; the absolute
        V^2/Hz matters to anyone integrating `S_v` to a power, and the
        scale is kept rather than doubled because it is a return value
        that callers may already divide by `|X_1|^2` themselves.  `S_v`
        was documented as the one-sided PSD of the output
        voltage; `L_dBc` is that normalised to the harmonic's own power,
        in dBc/Hz.

        ⚠ NO SWEEP AND NO PER-FREQUENCY SOLVE.  Once the PSS waveform's
        Fourier coefficients and the scalar `c` are known, "we have an
        analytical expression that gives us the spectrum at any frequency.
        The computation of the spectrum is not performed separately for
        every frequency of interest."  Which also means it never meets the
        near-carrier singularity that a swept small-signal computation
        would, and never meets the 1/f sweep-grid trap — there is no sweep
        to place a point on.

        ⚠⚠ SCOPE: A SOURCE BEHIND A SLOW NODE (A2, resolved 2026-09-08).
        The Lorentzian uses the DC PPV, so for a noise source that reaches
        the core through a slow path (RC leg, tau >> T) it holds only
        BELOW the source's corner `T/(2 pi tau)`; above it the true skirt
        is this one scaled by the PPV-harmonic-weighted filter
        `sum_k |G_k|^2 F_k(f) / sum_k |G_k|^2 F_k(0)` (G_k the PPV entry's
        Fourier coefficients at the source node, F_k the path's transfer at
        k f0 + f), which is 1/1000 at 0.1 f0 on a one-RC-leg fixture with
        an asymmetric core AND tank loss (both needed for G_0 != 0; an
        ideal tank inductor shorts DC).  `c` is still right (the filter
        removes only high-frequency content); `pnoise` computes the true
        value at any offset, and a Monte Carlo of `c` cannot see it.
        ✅ SINCE 2026-09-14 THIS METHOD DOES TOO, by default:
        `frequency_aware=True` replaces `c` by `c(f)` from the
        frequency-aware PPV (`frequency_aware_diffusion`), which matches
        pnoise's PM content to 0.4 % at 1e-3 and 1e-2 f0 on that fixture
        (the DC-PPV Lorentzian: 0.73x and 0.027x), and to <= 2 % above
        f_amp on an orbit with AM-to-PM coupling.  `frequency_aware=False`
        is the closed form, one `c` for every offset, and costs no solve.  Measured against `pnoise` to four digits
        through the corner (test ..._behind_a_slow_node_...).

        ⚠ AND IT IS THE ONLY ROUTE THAT IS VALID BELOW THE CORNER.  A
        small-signal analysis cannot produce `L(f)` there however well
        conditioned it is: the excess phase is a Wiener process, its
        spectrum has a singularity at the origin and no physical meaning,
        and the finite value `L` attains comes from the NONLINEAR
        phase-to-voltage map — which is what this closed form carries.
        Reporting `S_phi` near the carrier instead is the mistake that
        object invites.
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
        ## not represent.  Measured (2026-09-14, Monte Carlo, which agreed
        ## with pnoise to 1 % at every harmonic): van der Pol C=4 Q=8, 2 f0 +
        ## 10 f_amp -- this returned ~0 against 5.6e-6 V^2/Hz.  Refused, as
        ## `am_pm` refuses the same case.  ⚠ NOT caught, and recorded instead:
        ## away from the fundamental the model also misses where a line DOES
        ## exist (asymmetric orbit, 2 f0: 0.40 of the Monte Carlo) -- use
        ## `pnoise` away from the fundamental.
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

        ⚠⚠ WHY IT EXISTS (2026-09-14).  The Lorentzian from `c` uses the DC
        PPV at every offset: a noise current is assumed to move the phase
        instantly.  Wherever part of that response goes THROUGH a slow mode --
        the amplitude mode on an orbit with AM-to-PM coupling, or a slow node
        in the source's path -- it is filtered above that mode's corner, and
        the DC-PPV Lorentzian over-states.  Measured against pnoise's PM
        content (itself Monte-Carlo-confirmed on the first case):

            van der Pol C=4 Q=8, a=0.30     0.3 / 1 / 3 / 10 f_amp
              c(f)/c                         0.939 0.649 0.371 0.308
              pm / 4 S_v(DC PPV)             0.942 0.647 0.366 0.302
              pm / 4 S_v(c(f))               1.003 0.998 0.986 0.980
            A2 slow node (tau/T=100)         1e-3 / 1e-2 f0
              pm / 4 S_v(DC PPV)             0.729 0.027
              pm / 4 S_v(c(f))               1.004 1.004
            symmetric control (a=0)          c(f)/c within 1e-3

        ⚠ Stationary WHITE sources only, like `diffusion_constant`.  Cost:
        one bordered adjoint GMRES per offset (0.25-0.7 s on these fixtures);
        the solve can fail to converge (Lai's own warning about eq. 23), and
        then this raises rather than returning the DC value silently.
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
        adds to the total, so this method UNDER-reports.  Relayed measurement
        against a commercial simulator's total noise, as excess over the
        phase-only prediction:

            offset     lam2 = 0.90            lam2 = 0.99
                       (f_amp 26.7 kHz)       (f_amp 2.55 kHz)
            100 Hz     -0.00 dB               -0.01 dB
            1 kHz      -0.00 dB               -0.54 dB
            10 kHz     -0.50 dB               -2.90 dB
            100 kHz    -3.11 dB               -3.27 dB

        ⚠⚠ AND THE VALID REGION SHRINKS AS `1/Q`, which makes this section 0
        again rather than a detail.  With `f_amp = -ln(lam2)/(2 pi T)` and
        `Q = -1/ln(lam2)`,

            f_amp = f0 / (2 pi Q)

        -- verified both ways at 26671.9 / 2544.2 / 253.3 Hz for
        `lam2 = 0.90 / 0.99 / 0.999`.  So the better the oscillator, the
        narrower the band in which its phase-only spectrum is the whole
        answer; at `lam2 = 0.999` it has collapsed below ~253 Hz.

        ⚠ THIS IS THE OPPOSITE SIGN FROM THE ERROR `PSS.ppv` ALREADY WARNS
        ABOUT.  That one says the instantaneous phase equation misses slow
        nodes which FILTER device noise, so phase noise is OVER-estimated.
        This one is a second, independent mechanism in which the phase-only
        answer is UNDER-estimated.  Both are live and they are not the same
        effect.
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

        ⚠ THE REST OF THIS CLASS TAKES A DIRECTION VECTOR AND THIS PAIR
        TOOK AN INTEGER, which is not a style difference -- it meant
        `am_pm` and `carrier_phasor` could not express a DIFFERENTIAL
        output at all.  `pnoise`, `adjoint_transfer_row` and
        `adjoint_sideband_row` all accept `d`; these did `int(output)`.
        For an oscillator the output of interest is very often
        differential, and for the coordinate-invariance property an AM/PM
        split has to have (Kaertner 1990 section 3.2) a
        reference-independent observable is the whole point.

        An integer is still accepted, so callers that name a node keep
        working; an array is contracted against the full waveform with the
        reference row reinserted.
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

        ⚠ THE ABSOLUTE MAGNITUDE ON AN OSCILLATOR IS SMALL FOR A REASON
        (established 2026-09-08, the three-leg chain).  These are the p = 0
        band of `am_pm_noise`: a source at BASEBAND `freq` reaching the
        carrier sideband.  A baseband current moves the PHASE through the
        PPV's DC coefficient (Hajimiri-Lee's c_0), and a half-wave-symmetric
        orbit -- odd nonlinearity, `u(t + T/2) = -u(t)` -- has none, so on
        such a fixture the rows measure a symmetry zero (6e-9 .. 1e-13,
        proportional to 1/freq and to mu), the same zero the coloured
        up-conversion gate records for Gamma.  Breaking the symmetry
        (`_lc_osc(a)`) lifts |m_pm| at 1e-3 f0 from 1.2e-8 to 46.6 (a =
        0.05) and 231 (a = 0.25) -- linear in `a`.  The DIRECT rows (source
        at f0 + freq, sideband 0) agree with `pnoise` at every offset, and
        the split lands on the externally certified Lorentzian.  So do not
        read a small `am_pm` on a symmetric oscillator as a defect: it is
        the 1/f^3 up-conversion coefficient, and it is zero there.
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

        ⚠ THE AUTONOMOUS CAVEAT IS RETIRED (three-leg chain, 2026-09-08).  On a
        free-running oscillator this split sits on the SAME absolute scale as
        `pnoise` (identity to 1e-12 / 1e-16) and as the externally certified
        `oscillator_spectrum` (`S_pm = 4 S_v` at every offset: the PM content
        of the pair IS the Lorentzian, 2 S_v per sideband), with `S_am` rising
        from ~0 below the AM corner `f0/(2 pi Q_lambda)` to `S_pm` above it
        (⚠ this line said `f0/(4 pi Q)` until 2026-09-09 -- a factor of two
        the docs session caught against this very function: the ratio is an
        exact Lorentzian `u^2/(u_c^2 + u^2)` in `u = offset/f0` with
        `u_c = 1/(2 pi Q_lambda)`, `Q_lambda = -1/ln|lambda_2|`, half-power
        0.5007 there and 0.20 at the old corner, at Q = 8 and 16, 240 and
        480 points; the old formula OVERSTATED the AM content at every
        offset, 2.5x at its own corner) -- so the
        pair total is 4 S_v there and 8 S_v far out.  The "~1e-12 rows" were
        `am_pm`'s p = 0 band on a half-wave-symmetric fixture: a symmetry
        zero, see `am_pm`.  Oscillator magnitudes from this are trustworthy.

        ⚠ WHAT A MEASUREMENT MUST BE TO BE COMPARED WITH `S_pm` (2026-09-09,
        after a 16-seed Monte Carlo campaign and a deterministic forward
        tone route; roadmap "Item 3, answered"): `S_pm` is PM BY QUADRATURE
        OF THE FUNDAMENTAL'S SIDEBANDS.  A "phase" read by a one-period
        demodulation of the fundamental leaks the other harmonics'
        sidebands through its boxcar (sinc(pi(1 - r)) ~ 0.1 in amplitude
        for the second harmonic's, which is ~ the orbit's asymmetry), and
        a phase read from zero crossings converts EVERY harmonic's
        sidebands; both drifted 6 % against this quantity as the asymmetry
        of a harmonic-rich orbit was swept while the forward LPTV response
        agreed with it to 1 %.  So compare `S_pm` with the fundamental's
        sideband PM (a spectrum analyser's sidebands around f0, or the
        forward-tone gate `test_pnoise_oscillator_pm_matches_a_forward_tone_
        transient_with_no_adjoint`), never with a demodulated or
        crossing-time phase, and BAND WITH BAND: a source behind a slow node
        has an in-band spectrum that is not 1/r^2 (its slow/core ratio
        swings 1.16 -> 0.87 across 0.08-0.15 f0), so a band mean and a
        point value differ by ~4 % there.
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
        ## ⚠⚠ THE SPLIT IS TAKEN IN THE CARRIER'S FRAME, NOT THE TIME ORIGIN'S
        ## (defect reported by a peer session and reproduced 2026-09-14).  AM
        ## is the envelope component ALONG the carrier phasor, so `a + conj(b)`
        ## is right only for a cosine-phased carrier.  Until this rotation the
        ## answer depended on where t = 0 sat: a driven diode gave am/pm =
        ## 0.305 / 3.28 / 0.651 at drive phase 0 / 90 / 37 degrees, and
        ## 3.316 at all three once rotated, `S_am + S_pm` unchanged to 1e-15
        ## (the identity cannot see it -- `|a_r|`, `|b_r|` are `|a|`, `|b|`).
        ## `am_pm` never had it: it divides by the COMPLEX carrier phasor.
        ## With no carrier at this harmonic the phase is undefined and the
        ## split is left unrotated, as `am_pm` refuses the same case.
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

        MEASURED on van der Pol, offsets from 0.3 down to 1e-9 of `f0`:

            offset/f0    sigma_min(plain)   sigma_min(bordered)
            3e-01           5.68e-01            1.17e-01
            1e-03           2.61e-03            2.04e-01
            1e-06           2.61e-06            2.04e-01
            1e-09           2.61e-09            2.04e-01

        The plain operator tracks the offset over nine decades; the
        bordered one is FLAT.  The two solutions agree to 5.7e-12 where the
        plain solve is still trustworthy, and their disagreement grows as
        `1/df` -- that is the PLAIN solve losing digits, not this one.

        ⚠ IT STILL DIVERGES AT AN EXACT HARMONIC, and it should: `1/(1 -
        alpha)` is then a division by zero, and the physical response is
        unbounded.  What changes is that every offset NEAR a harmonic is
        now well conditioned, which is where phase noise is measured.

        `transposed` solves `(I - alpha M^T) x = b`, whose null space is
        spanned by `v` and whose left null space is spanned by `u`, so the
        borders swap.
        """
        import scipy.sparse.linalg as spla
        fp = pss.factored_period()
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
        ## ⚠ THE BORDER VECTORS ARE NORMALISED (2026-09-22): the recovery
        ## `y = w + s col / (1 - alpha)` is scale-free in exact arithmetic,
        ## but GMRES sees the bordered MATRIX, and with the tangent in V/s
        ## (1e6 on a comparator oscillator) against the PPV in s/V (1e-6)
        ## its condition number was 1.8e9 where `I - alpha M` alone was
        ## 2.4 -- the forward solve returned with a 5.5e-4 residual and
        ## the transposed one did not converge at all.  Unit vectors put
        ## the bordering at the operator's own conditioning; `s` absorbs
        ## the scale.
        col = col / max(float(np.linalg.norm(col)), 1e-300)
        row = row / max(float(np.linalg.norm(row)), 1e-300)
        mv = (fp.matvec_transposed if transposed else fp.matvec)
        ## ⚠ ON A STAGED SOLVE THE POLE IS THE TOTAL MAP'S (2026-09-22,
        ## events phase B): `u`, `v` are the null vectors of `I - M_tot`,
        ## `M_tot = M + P_theta dtheta/dx_0`, and the fixed-grid `M` has no
        ## unit multiplier at all (|lambda - 1| = 0.99 on the comparator
        ## oscillator) -- bordered with the total map's vectors it read
        ## the sideband response 0.3-400x off the exact one.  The
        ## operator here is the total map's.
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
        ## ⚠ REFINED ON THE PLAIN OPERATOR WHERE THAT IS WELL CONDITIONED
        ## (2026-09-22): the recovery assumes `u`, `v` are EXACT null
        ## vectors of `I - M`.  On a staged solve the discrete total map's
        ## unit multiplier is displaced by O(h) (8e-4 at 200 points on the
        ## comparator oscillator, 4e-4 at 400; an unstaged radau map sits
        ## at 1e-11), and the recovered `y` then misses the true operator
        ## by that much over `|1 - alpha|`: 5.5e-4 at 0.3 f0, 0.14 at
        ## 1.001 f0.  The plain operator is well conditioned there (2.4 at
        ## 0.3 f0, 310 at 1.001 f0, 2e3 at 1.00001 f0 -- its pole sits
        ## where the DISCRETE multiplier is, not at alpha = 1), so the
        ## deflated answer seeds a plain correction on its own residual,
        ## kept only if it lowers the residual; below
        ## `DEFLATION_REFINE_MIN` the deflated answer stands.  ⚠ This makes
        ## the forward and adjoint solves the discrete operator's own,
        ## dual-consistent -- and near a harmonic the discrete operator's
        ## answer carries the multiplier's displacement: 15.7 % off the
        ## exact forced response at 1.001 f0 on the 200-point staged
        ## comparator oscillator (0.2 % at 0.3 f0 and 1.7 f0), where the
        ## unrefined recovery, which carries the pole analytically at
        ## alpha = 1, happened to read 0.2 % but is not dual-consistent
        ## (4.5e-4 at 0.3 f0) and has its own O(h/|1 - alpha|) split
        ## error.  The item is the staged map's multiplier (roadmap E8).
        ## On an unstaged oscillator the residual is already at the
        ## tolerance and nothing happens.
        if abs(denom) >= self.DEFLATION_REFINE_MIN:
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
