"""`PAC`, the periodic small-signal and noise analyses over a `PSS` operating
point, and `SidebandResponse`.

`PAC` keeps its small-signal core here -- the forced and adjoint sideband
solves, the mixer response, AM/PM and the deflated solve -- and inherits its
noise themes from the `_pac_*` modules:

    _pac_sources.py    the noise sources: CY on the orbit, the colour models,
                       their fits and roots, the guards, the quadrature
    _pac_pnoise.py     pnoise and its folds, am_pm_noise, band_spread
    _pac_lyapunov.py   the periodic covariance, event jitter
    _pac_osccov.py     an oscillator's covariance and edge jitter
    _pac_sampled.py    sampled noise, its variance, jitter metrics
    _pac_modal.py      the Floquet-mode (orbital, modal) spectra
    _pac_phase.py      phase noise: diffusion, phase_psd, the lineshape

History: `doc/shooting_history.md`, `pac` (module level).
"""
import numpy as np
import weakref
from pycircuit.circuit.analysis import Analysis
from pycircuit.circuit.analysis import Parameter
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.circuit import gnd
import pycircuit.circuit.analysis as analysis
from ._factored import dense_map
from ._numerics import _arnoldi_gmres
from ._numerics import _output_weights, output_index
from ._numerics import sweep_frequency, sweep_offset
from ._numerics import freq_analysis
from ._pac_lyapunov import _LyapunovCovariance
from ._pac_modal import _ModalSpectra
from ._pac_osccov import _OscillatorCovariance
from ._pac_phase import _PhaseNoise
from ._pac_pnoise import _DrivenNoise
from ._pac_sampled import _SampledNoise
from ._pac_sources import _NoiseSources
from .events import EventColumns
from pycircuit.circuit.simwarnings import AccuracyWarning, warn


class _SidebandFamily:
    """The adjoint sideband rows sharing one OUTPUT frequency
    (`PAC._sideband_family`): `parts` are `(times, couplings, scale)` from
    its reverse passes, and the row of sideband `l` at input frequency
    `fin` is ``-sum scale * exp(2j pi fin times) @ couplings``."""
    __slots__ = ('_N', '_T', '_m', '_pac', '_parts', '_pss', 'f_out',
                 'matvecs')

    def __init__(self, pac, pss, f_out, T, N, m, parts, matvecs):
        self._pac, self._pss, self.f_out = pac, pss, f_out
        self._T, self._N, self._m = T, N, m
        self._parts, self.matvecs = parts, matvecs

    def row(self, l, fin=None):
        """The row of sideband `l`: every source at `fin` (default
        `f_out - l f0`) to the output's sideband `l`."""
        fin = (self.f_out - float(l) / self._T) if fin is None else float(fin)
        if abs(int(l)) > self._N // 2:
            raise ValueError(
                f'PAC: sideband {int(l)} is above the grid\'s Nyquist (|l| <= '
                f'{self._N // 2} at {self._N} points per period). Nothing can '
                'alias down from above the maximum frequency the grid '
                'represents, so this is not a tolerance to relax -- use a '
                'finer period grid.')
        self._pac._check_harmonic(self._pss, fin, 'the sideband row')
        out = np.zeros(self._m, dtype=complex)
        for t, c, sc in self._parts:
            if len(t):
                out = out - sc * (np.exp(2j * np.pi * fin * t) @ c)
        return out


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


class PAC(_NoiseSources, _DrivenNoise, _LyapunovCovariance,
          _OscillatorCovariance, _SampledNoise, _ModalSpectra,
          _PhaseNoise, Analysis):
    """Small-signal analysis over a periodic operating point, matrix-free.

    The noise, jitter and phase-noise methods' conventions -- which are
    one-sided, which 0.5x one-sided or two-sided, which frequency argument
    is an offset, what each returns -- are tabulated in
    `doc/pac_noise_conventions.md`.

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

    def solve(self, pss, freqs, refnode=gnd, recycle=True, sweeptype=None,
              relharmnum=None):
        """Sideband response at each frequency in `freqs`.

        `freqs` are the SOURCE's frequencies, by a commercial RF simulator's
        sweep rule: `sweeptype='absolute'` the frequencies themselves,
        `'relative'` offsets `relharmnum * f0 + freqs` (`relharmnum` default
        1), `None` relative on an AUTONOMOUS PSS and absolute on a driven
        one (`_numerics.sweep_kind`).  ⚠ Until 2026-09-29 absolute on every
        PSS.  The result's sweep values are the absolute OUTPUT frequencies.

        `pss` must be a CONVERGED `PSS` -- the periodic operating point is
        what this linearises about, and there is no meaningful small-signal
        answer about a non-solution.  `PSS.factored_period()` enforces it.

        `recycle` shares one Krylov subspace across the sweep, which is
        where the sweep's cost goes; see `_solve_subspace`.

        Returns a `CircuitResult` (also `self.result`) swept over the OUTPUT
        frequency: every sideband ``f + k f0`` of every sweep point, merged
        and sorted, a negative one folded to ``|f + k f0|`` with its
        coefficient conjugated.  Also sets `self.time_response`, per sweep
        point ``(times, y)`` the complex response at the grid's nodes (a
        time-domain reading comes from here, not from summing sidebands),
        and `self.event_shifts`, per sweep point the landed crossings'
        small-signal shifts (None without state events).

        History: `doc/shooting_history.md`, `PAC.solve`.
        """
        toolkit = self.toolkit
        ## ⚠ ONE HOST: the map, its event columns, its period and orbit all
        ## the monodromy twin's (a trap/euler oscillator's radau twin; `self`
        ## otherwise).  Until 2026-09-30 the map was the twin's and the event
        ## columns the run's -- a trap staged oscillator crashed on the two
        ## grids (270 against 236 nodes) -- and a relative sweep offset from
        ## the run's carrier.
        pss = pss.monodromy_twin()
        freqs = np.atleast_1d(np.asarray(
            sweep_frequency(pss, np.asarray(freqs, dtype=float), sweeptype,
                            relharmnum, 'solve'), dtype=float))
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
            warn(
                'PAC: this operating point was solved on the PLAIN path '
                'with a manufacturing step (method=%r, x0_unknown=False). '
                'The manufacturing step carries no small-signal source, so '
                'the response is FIRST order in the timestep whatever the '
                "method's own order -- measured 2.00x per doubling against "
                '4.00x for the same run with x0_unknown=True. The answer is '
                'not wrong, it is one order less accurate than the '
                'trajectory it came from. Re-solve with x0_unknown=True, or '
                "with method='gear', to get the method's own order."
                % pss.par.method, AccuracyWarning)

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
        fixed-time correction of the node responses -- `PAC.solve`'s core,
        shared with `_coloured_covariance`.  Sets `deflated` and `matvecs`.

        History: `doc/shooting_history.md`, `PAC._forced_responses`."""
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
        _ev = EventColumns.of(pss)
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

        `recycle_tol` is the GMRES tolerance of the transposed solve
        (default `KRYLOV_FACTOR * reltol`), and binds only above
        `FLOQUET_DENSE_LIMIT`: below it the solve is direct.

        History: `doc/shooting_history.md`, `PAC.adjoint_transfer_row`.
        """
        output = output_index(pss, output)
        import scipy.sparse.linalg as spla
        ## (one host, as `PAC.solve`: the twin's map WITH the twin's event
        ## columns and period -- the run's were read beside the twin's map
        ## until 2026-09-30)
        pss = pss.monodromy_twin()
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
        tol = max(tol, 1e-14)
        A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
        ## the same pole as in `solve` and `adjoint_sideband_row`: deflated
        ## on an oscillator, plain (and cheaper) on a driven circuit
        _autonomous = bool(getattr(pss, 'autonomous', False))
        self.deflated = _autonomous
        _ev = EventColumns.of(pss)
        if _ev is None:
            xa = self._pole_solve(pss, A, alpha, d, tol, 'the adjoint solve')
            self.matvecs = count[0]
            return alpha * pss._forced_replay_transposed(fp, freq, xa)
        ## ⚠ ON A STAGED SOLVE THE ROW IS BORDERED, as `adjoint_sideband_row`'s
        ## (`_sideband_family`): the transpose of `PAC.solve`'s bordered
        ## system with the output `d . y_0` -- read at node 0, whose time the
        ## events do not move (`g_theta` over the fixed-time columns, zero
        ## there) -- and the event rows' term as a second reverse pass.
        ## Unbordered until 2026-10-01 (the review's D3), it missed the
        ## bordered forward solve by 2.8e-4 (radau) / 1.6e-3 (trap) / 3.9e-2
        ## (gear) on the PWM loop at 60 points.
        N = len(fp.steps)
        cn = np.zeros(N, dtype=complex)
        cn[0] = 1.0
        _Pkf = self._fixed_time_event_columns(pss)[0]
        g_theta = EventColumns.g_theta(cn, _Pkf, d[:m], N)
        if _autonomous:
            z = self._deflated_solve(
                pss, alpha, np.asarray(d, dtype=complex)
                + np.asarray(_ev.dth, dtype=float).T @ g_theta,
                transposed=True, tol=tol)
            zeta = _ev.collapsed_zeta(g_theta, alpha, z)
        else:
            def _solve_adj(b, k):
                return self._adjoint_solve(
                    pss, fp, alpha, b, A, tol,
                    'the adjoint solve' if k is None
                    else f'the bordered adjoint solve, event {k}',
                    total=False)
            z, zeta = _ev.bordered_adjoint(_solve_adj, d, g_theta, alpha)
        self.matvecs = count[0]
        row = alpha * pss._forced_replay_transposed(fp, freq, z)
        t_e, c_e, _g = pss._reverse_points(fp, freq,
                                           extra=_ev.injection_dict(zeta))
        if len(t_e):
            row = row - np.exp(2j * np.pi * float(freq) * t_e) @ c_e
        return row

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
        return self._adjoint_solve(pss, pss._state_map(), alpha, rhs, A,
                                   tol, what, total=False)

    def _dense_period(self, pss, fp, n, total=True):
        """The period map as a dense matrix -- the TOTAL one on a staged
        solve (`EventColumns.total_matrix`) unless `total=False` -- or None
        above `FLOQUET_DENSE_LIMIT`, where the solves stay matrix-free."""
        if n > pss.FLOQUET_DENSE_LIMIT:
            return None
        M = dense_map(fp, n)
        _ev = EventColumns.of(pss, n) if total else None
        return M if _ev is None else _ev.total_matrix(M)

    def _adjoint_solve(self, pss, fp, alpha, b, A, tol, what, total=True):
        """``(I - alpha M^T) x = b`` on a DRIVEN circuit (no pole): directly
        from the dense map below `FLOQUET_DENSE_LIMIT` (2026-09-30; exact to
        rounding), GMRES on `A` above it."""
        b = np.asarray(b, dtype=complex).ravel()
        Md = self._dense_period(pss, fp, len(b), total=total)
        if Md is None:
            return self._gmres_checked(A, b, tol, what)
        return np.linalg.solve(np.eye(len(b)) - alpha * Md.T, b)

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
        output = output_index(pss, output)
        ## (one host, as `PAC.solve`: the twin's map WITH the twin's event
        ## columns and period -- the run's were read beside the twin's map
        ## until 2026-09-30)
        pss = pss.monodromy_twin()
        fp = pss._state_map()

        self._check_circuit(pss)
        self._check_harmonic(pss, freq, 'the sideband row')
        m = pss.cir.n - 1
        T = float(fp.T)
        N = len(fp.steps)
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

        rows = np.zeros((len(ls), m), dtype=complex)
        nmv = 0
        for li, l in enumerate(ls):
            ## (one family per OUTPUT frequency `l f0 + freq`; `pnoise`,
            ## `am_pm_noise` and `mixer_response` share one family across
            ## every sideband of theirs -- `_sideband_family`)
            fam = self._sideband_family(pss, float(l) / T + float(freq), output)
            rows[li] = fam.row(int(l), float(freq))
            nmv += fam.matvecs
        self.matvecs = nmv
        return rows

    def _sideband_family(self, pss, f_out, output):
        """Every adjoint sideband row with OUTPUT frequency `f_out` -- the row
        of sideband `l` has input frequency `f_out - l f0` -- as a
        `_SidebandFamily`, whose `row(l, fin)` is a phase-weighted sum.

        ⚠ ONE PASS, ONE SOLVE, ONE PASS FOR THE WHOLE FAMILY.  The output
        functional is injected with phase ``exp(-j (l w0 + w_in) t)`` =
        ``exp(-2j pi f_out t)`` and ``alpha = exp(-2j pi f_in T)`` =
        ``exp(-2j pi f_out T)`` -- neither depends on `l`, so the reverse
        pass, its final costate `g`, the pole/bordered solve for `z` and
        the event rows are the family's; only the SOURCE's phase at each
        injection point is the row's, and each pass collects its couplings
        (`_reverse_points`) instead of applying them at one frequency.
        Until 2026-09-30 `pnoise` repeated all of it per sideband (83 % of
        its time); the costate agreed across `l` to 9.4e-15.

        History: `doc/shooting_history.md`, `PAC.adjoint_sideband_row`."""
        import scipy.sparse.linalg as spla
        output = output_index(pss, output)
        pss = pss.monodromy_twin()
        fp = pss._state_map()
        self._check_circuit(pss)
        m = pss.cir.n - 1
        n = fp.width
        T = float(fp.T)
        tms = np.asarray(fp.times, dtype=float)
        N = len(fp.steps)
        alpha = np.exp(-2j * np.pi * float(f_out) * T)
        d = _output_weights(output, m)
        count = [0]

        def _mv(v):
            count[0] += 1
            return np.asarray(v) - alpha * fp.matvec_transposed(v)

        A = spla.LinearOperator((n, n), matvec=_mv, dtype=complex)
        tol = max(self.KRYLOV_FACTOR * pss.par.reltol, 1e-14)
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
        ## every node -- one fold for every kind (`_reverse_points`,
        ## verified vs forward driven solves and the bespoke trbdf2 fold)
        t_f, c_f, g = pss._reverse_points(fp, f_out, d)
        parts = [(t_f, c_f, 1.0)]
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
        ## solve is 10-15 % off (driven) and the row misses the forward
        ## solve by 8 % on gear's staged oscillator.
        ## History: `doc/shooting_history.md`, `PAC.adjoint_sideband_row`.
        ## ⚠ EVERY MAP KIND: the plain map (trap/euler) stayed unbordered
        ## until 2026-09-30 -- the kinds were admitted one by one as each
        ## was verified and the plain one never was -- and its row missed
        ## the bordered forward solve by 1.9e-3 (trap) / 9e-3 (euler).
        _autonomous = getattr(pss, 'autonomous', False)
        _ev = EventColumns.of(pss)
        if _ev is not None:
            _wq = pss._period_quadrature(fp)
            cn = np.array([np.exp(-2j * np.pi * float(f_out) * tms[j])
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
                    return self._adjoint_solve(
                        pss, fp, alpha, b, A, tol,
                        ('the adjoint solve at %.6g Hz' % f_out) if k is None
                        else ('the bordered adjoint solve, event %d' % k),
                        total=False)
                z, zeta = _ev.bordered_adjoint(_solve_adj, g, g_theta, alpha)
            t_e, c_e, _g2 = pss._reverse_points(
                fp, f_out, extra=_ev.injection_dict(zeta))
            parts.append((t_e, c_e, 1.0))
        else:
            z = self._pole_solve(pss, A, alpha, g, tol,
                                 'the adjoint solve at %.6g Hz' % f_out)
        t_z, c_z, _g3 = pss._reverse_points(fp, f_out, lam0=z)
        parts.append((t_z, c_z, alpha))
        return _SidebandFamily(self, pss, float(f_out), T, N, m, parts,
                               count[0])

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
        output = output_index(pss, output)
        self._check_circuit(pss)
        f0 = 1.0 / float(pss.period)
        ls = [int(l) for l in sidebands]
        ins, rows = [], []
        fam = None
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
            ## (every band lands on `f_out`: one family, `_sideband_family`)
            if fam is None:
                fam = self._sideband_family(pss, float(f_out), output)
            row = fam.row(int(l), f_in)
            ins.append(f_in)
            rows.append(np.asarray(row).reshape(-1))
        return SidebandResponse(f_out, f0, ls, ins, np.array(rows))

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
                    ## (at the solve's `epar`, as `_k_at` reads `u`: a source
                    ## may depend on the temperature -- the review's D6)
                    ud = np.delete(np.asarray(pss.cir.dudt(t, epar=pss.epar,
                                                           analysis=analysis),
                                              dtype=float).ravel(), pss.irefnode)
                    A[alg, :] = G[alg, :]
                    b[alg] = -ud[alg]
                out[j] = np.linalg.solve(A, b)
            except (np.linalg.LinAlgError, ValueError):
                ok = False
                break
        if ok:
            return out
        warn(
            'PAC._orbit_rate: the DAE derivative could not be assembled at a '
            'node (a singular differential/algebraic split -- an index above '
            'one?); falling back to the three-node stencil, which is second '
            'order in the step and one-sided at a landed event.', AccuracyWarning)
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
        ev = EventColumns.of(pss)
        th = np.asarray(pss._state_event_fracs, dtype=float)
        _fr, hsens, _nd = pss._event_remap(
            np.asarray(pss._grid_fracs, dtype=float), th, th, float(pss.period))
        tau = np.vstack((np.zeros((1, len(th))),
                         np.cumsum(np.asarray(hsens, dtype=float), axis=0)))
        xdot = self._orbit_rate(pss, ev['nodes'])
        Pk = np.asarray(ev['Pk_nodes'], dtype=float)
        return Pk - xdot[:, :, None] * tau[:, None, :], tau, xdot

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
        order of `_stage_times`: each step's own, as the factored walk
        solved them (`_StageStep.Ys`; a GLM's steps list their startups'
        substages too).
        History: `doc/shooting_history.md`, `PAC._stage_states`."""
        if fp.is_glm:
            return [y for st in fp.step_objects()
                    for y in st.injection_states()]
        return [y for st in fp.steps for y in st.Ys]

    def _stage_pass(self, pss, fp, lam0, seed=None):
        """One reverse pass of a stage period map (`dirk` or `full`, or a
        GLM's map on the state): returns the final costate and the coupling
        vectors, one per injection point (`_stage_times`; `h sum_i A_ik p_i`
        per stage of a stage method) -- the sensitivity of the costate's
        functional to a unit source there is minus that (see
        `_reverse_points`, whose loop this is).  `seed = (k0, v)` adds `v`
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

    def carrier_phasor(self, pss, output, harmonic=1):
        """The `harmonic`-th Fourier coefficient of the steady-state output.

        Computed here rather than taken from `fpss`, whose spectrum is RMS
        and energy-folded -- correct for reporting a magnitude and useless
        for a phasor, since folding discards the phase the AM/PM split is
        made of.
        """
        output = output_index(pss, output)
        times, _X = pss.waveform
        row = self._output_waveform_row(pss, output)
        t = np.asarray(times, dtype=float)[:-1]
        v = row[:len(t)]
        w0 = 2.0 * np.pi / float(pss.period)
        ## ⚠ a Fourier INTEGRAL: `1/N` is its quadrature only on a uniform
        ## grid (measured 7.5 % off and not converging on a 3:1 one).  ⚠ The
        ## WAVEFORM's own grid: `factored_period()` is the twin's on a trap
        ## or euler oscillator, and its weights against the run's samples
        ## put the carrier 90 % off on a staged trap one (until 2026-10-01,
        ## the review's D7)
        _wq = pss._times_quadrature(np.asarray(times, dtype=float))
        if _wq is None or len(_wq) != len(t):
            return complex(np.sum(v * np.exp(-1j * harmonic * w0 * t)) / len(t))
        return complex(np.sum(v * np.exp(-1j * harmonic * w0 * t) * _wq))

    def am_pm(self, pss, freq, output, harmonic=1, sweeptype=None):
        """AM and PM modulation indices at `harmonic`, per noise/signal source.

        Returns `(m_am, m_pm)`, each a row of length `m`: the modulation a
        unit source at reduced coordinate `i`, driven at `freq`, imposes on
        the `harmonic`-th harmonic of the output.

        `freq` by a commercial RF simulator's sweep rule with `harmonic` as
        the reference, as `am_pm_noise`: `sweeptype='relative'` the offset
        from `harmonic*f0` (the source's own, baseband, frequency),
        `'absolute'` the upper output sideband's frequency, `None` relative
        on an AUTONOMOUS PSS and absolute on a driven one.  ⚠ Until
        2026-09-29 an offset on every PSS.  Below, `freq` is the offset.

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
        output = output_index(pss, output)
        ## (one host, as `am_pm_noise`: the carrier, its phasor and the rows
        ## all the map's the rows are solved on)
        pss = pss.monodromy_twin()
        freq = sweep_offset(pss, freq, sweeptype, harmonic, 'am_pm')
        C = self.carrier_phasor(pss, output, harmonic)
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
                % (harmonic, abs(C), scale))
        upper = self.adjoint_sideband_row(pss, freq, output, harmonic)[0]
        lower = self.adjoint_sideband_row(pss, -freq, output, harmonic)[0]
        return self.am_pm_indices(upper / C, lower / C)

    ## `|1 - alpha|` below which the deflated answer is left as recovered:
    ## nearer a harmonic than this the plain operator is singular at the
    ## arithmetic (an unstaged radau map's multiplier sits at 1e-11 from 1)
    DEFLATION_REFINE_MIN = 1e-8

    def _deflated_solve(self, pss, alpha, b, transposed=False, tol=None,
                        parts=False):
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
        ## ⚠ THE BORDER, ONCE PER SOLVED ORBIT: `v` and `u` depend only on
        ## the converged map, and `ppv()` costs 0.22 s per call.  Kept per
        ## (pss, state map, the `ppv` that made them); a re-solve builds a
        ## new map and misses.  ⚠ THE PRODUCER IS PART OF THE KEY: a caller
        ## who replaces `pss.ppv` (the border-sensitivity test injects
        ## perturbed vectors that way) must be consulted.
        ## History: `doc/shooting_history.md`, `PAC._deflated_solve`.
        _ppv_fn = getattr(pss.ppv, '__func__', pss.ppv)
        _border = getattr(self, '_deflation_border', None)
        ## (the analysis held WEAKLY: a PAC kept the last PSS alive after
        ## its caller dropped it -- the review's M2)
        if (_border is not None and _border[0]() is pss and _border[1] is fp
                and _border[4] is _ppv_fn):
            v, u = _border[2], _border[3]
        else:
            _v, info = pss.ppv()
            v = np.asarray(_v, dtype=float)
            u = np.asarray(info['tangent_pair'], dtype=float)
            _pr = ((lambda _p=pss: _p) if isinstance(pss, weakref.ProxyTypes)
                   else weakref.ref(pss))
            self._deflation_border = (_pr, fp, v, u, _ppv_fn)
        vu = float(v @ u)
        ## (relative to the vectors' sizes: an absolute 1e-300 read nothing
        ## -- the review's X5)
        if abs(vu) <= 1e-14 * float(np.linalg.norm(v) * np.linalg.norm(u)):
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
        Md = self._dense_period(pss, fp, n)
        if Md is not None:
            ## ⚠ DENSE BELOW THE LIMIT: the bordered system from the cached
            ## map (`dense_map`), solved directly -- 0.14 ms against 0.55 s
            ## of GMRES at the review's measurement (2026-09-30), and exact
            ## to rounding rather than to the Krylov tolerance
            Bd = np.zeros((n + 1, n + 1), dtype=complex)
            Bd[:n, :n] = np.eye(n) - alpha * (Md.T if transposed else Md)
            Bd[:n, n] = col
            Bd[n, :n] = row
            z = np.linalg.solve(Bd, rhs)
        else:
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
        ## `parts`: also the BOUNDED part, `w` -- the response with the pole's
        ## `s col / (1 - alpha)` never added, which `_transverse_responses`
        ## replays so the 1/offset term cannot enter its projection
        wb = w.copy() if parts else None
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
        ## ⚠⚠ NOT ON A NORDSIECK GLM'S MAP.  Its map opens with the startup,
        ## which breaks the discrete phase symmetry: the unit multiplier sits
        ## ``eta = O(h^p)`` off 1.  Refined, the answer is the discrete
        ## operator's and misses the physical one by ``eta / (2 pi r)``, `r`
        ## the offset in units of f0.  Unrefined -- the pole carried
        ## analytically -- it is O(h^p) at every offset, and forward and
        ## adjoint agree to O(eta) rather than the arithmetic.
        ## History: `doc/shooting_history.md`, `PAC._deflated_solve`.
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
                    if parts:
                        ## the correction's own bounded part (small: no
                        ## cancellation against the pole)
                        wb = wb + dy - col * (complex(row @ dy)
                                              / complex(row @ col))
        return (y, wb) if parts else y

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
