"""`ProbeShooting`: the probe (harmonic-balance-like) shooting of B5.
"""
import numpy as np
from pycircuit.circuit.analysis import remove_row_col
from pycircuit.circuit.circuit import gnd
from .pac import PAC
from .pss import PSS


class ProbeShooting:
    """Bizzarri's probe-based shooting -- oscillator amplitude and frequency
    from a DRIVEN solve, plus the 2x2 power-flow instability screen.

    A periodic voltage source of amplitude ``A`` and frequency ``f`` is placed
    across ``node``, and ``(A, f)`` is solved so that the probe's OWN CURRENT
    vanishes.  At that point the probe sources nothing and can be removed
    without changing the steady state, so it is the oscillator's own solution.

    The probe makes the circuit NON-AUTONOMOUS: the period is known (``1/f``),
    so there is no phase condition, no free-period unknown, and no ``T = 0``
    trivial root for a seed below the fundamental to fall into.

    ⚠ IT IS NOT A CONVERGENCE AID, AND THE PAPER SAYS SO ITSELF.  On its own
    flagship high-Q Pierce example the authors report *"it is easy to assign a
    tentative current to the Ls inductor ... and obtain convergence in a few
    iterations (we did this with conventional SH)"*.  What this buys is the
    SWEEP: unstable limit cycles, coexisting solutions, and a stability screen.

    ⚠⚠ ONE TONE GIVES A DESCRIBING-FUNCTION SOLVE, NOT THE ORBIT.  A single
    tone forces a SINUSOID, so a non-sinusoidal orbit can only null the probe's
    FUNDAMENTAL current.  Measured on van der Pol, the frequency error is
    QUADRATIC in harmonic content -- ``df/f = 4.0 THD^2`` to 3% across two
    decades -- and the probe sits at the LC resonance at every ``mu`` because
    ``mu(u - u^3/3)`` is odd and memoryless, so its describing function shifts
    no phase.  ``harmonics=K`` forces K tones and nulls K harmonics, which is
    harmonic balance with a shooting inner solve: at K=3 the van der Pol error
    falls from 5.9e-02 to 3.8e-04.

    ⚠⚠ PROBE PLACEMENT IS CIRCUIT-SPECIFIC, AND ITS FAILURE IS NOT A SOLVER
    FAILURE.  Forcing a node fixes every state the source reaches; any state it
    does NOT reach whose DC level is then unconstrained makes the shooting
    Jacobian SINGULAR, because a whole family satisfies periodicity.  Measured
    across van der Pol's only node with no series resistance: periodicity error
    **2.11e-15** -- already a periodic solution -- reported as
    ``converged = False``.  :meth:`degenerate_placement` names that pairing.
    """

    def __init__(self, factory, node, refnode=gnd, method='gear',
                 reltol=1e-10, npts=300, maxiterations=30, phase=90.0,
                 harmonics=1, tones=None, warm_start=True):
        """`factory()` must return a FRESH circuit WITHOUT the probe.

        `tones` overrides `harmonics` with an explicit list of harmonic
        numbers (``[1, 3, 5]`` to skip the even ones).  ⚠ Do not choose it by
        eye -- see :meth:`even_harmonic_content`.
        """
        self.factory = factory
        self.node = node
        self.refnode = refnode
        self.method = method
        self.reltol = float(reltol)
        self.npts = int(npts)
        self.maxiterations = int(maxiterations)
        self.phase = float(phase)
        self.tones = ([int(k) for k in tones] if tones is not None
                      else list(range(1, int(harmonics) + 1)))
        if self.tones[0] != 1:
            raise ValueError('the fundamental (tone 1) must be present, got %r'
                             % (self.tones,))
        self.warm_start = bool(warm_start)
        self._x0 = None
        self.evaluations = 0

    ## ---- circuit construction -------------------------------------------

    def _build(self, f, amps, phases):
        """The oscillator with K probe sources IN SERIES across `node`.

        In series they share ONE physical current (KCL ties them at the
        intermediate nodes), which is the quantity to null -- so K tones need
        no new element type, only K sources and K-1 nodes.
        """
        from pycircuit.circuit.elements import VSin
        cir = self.factory()
        prev = self.node
        n = len(self.tones)
        for j, k in enumerate(self.tones):
            nxt = self.refnode if j == n - 1 else cir.add_node('__probe_n%d' % j)
            ## ⚠⚠ `vac=0` EXPLICITLY.  `VS.vac` DEFAULTS TO 1, not 0, so every
            ## probe source in the chain would be AC-excited at once -- and in
            ## series they share one current, so the PAC response came back
            ## exactly K times too large (measured ratios 1, 2, 3 at K = 1, 2,
            ## 3).  K=1 validated because there was nothing to contaminate it,
            ## which is why the diagonal alone could not catch this.
            cir['__probe%d' % j] = VSin(prev, nxt, va=float(amps[j]),
                                        freq=float(k) * float(f),
                                        phase=float(phases[j]), vac=0.0)
            prev = nxt
        return cir

    def _probe_row(self, cir):
        """The global row of the probe chain's branch current."""
        rows = cir.elementnodemap['__probe0']
        el = cir['__probe0']
        if len(el.branches) != 1:
            raise ValueError(
                'ProbeShooting expects the probe to carry exactly one branch, '
                'this one carries %d.' % len(el.branches))
        return int(rows[-1])

    ## ---- evaluation ------------------------------------------------------

    def _spectrum(self, f, amps, phases, upto=None):
        """`(I, pss)` -- the probe current's harmonic phasors.

        `upto` defaults to the configured tones; pass a larger count to look at
        harmonics the solve is NOT nulling, which is what
        :meth:`even_harmonic_content` needs.
        """
        import warnings as _w
        cir = self._build(f, amps, phases)
        row = self._probe_row(cir)
        pss = PSS(cir, method=self.method, reltol=self.reltol)
        T = 1.0 / float(f)
        ## ⚠ WARM START: the finite-difference columns perturb a parameter by
        ## ~1e-5, so the trajectory barely moves and solving each from cold is
        ## waste.  Measured 2.09 s cold against 1.31 s warm -- 1.60x -- with
        ## the answers agreeing to 3e-15.  A pure accelerator: it changes which
        ## iterate the Newton starts from and nothing else, and carries no
        ## assumption about the circuit.
        ## ⚠ `PSS.solve` takes the REDUCED state (length n-1), not the full
        ## one -- passing `X[:, 0]` verbatim makes a solved-history run build a
        ## 2(n-1) pair against an n-length vector and raise on the shapes.
        x0 = self._x0 if (self.warm_start and self._x0 is not None
                          and len(self._x0) == cir.n - 1) else None
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / self.npts, x0=x0,
                      maxiterations=self.maxiterations)
        self.evaluations += 1
        X = np.asarray(pss.waveform[1], dtype=float)
        if pss.converged:
            self._x0 = np.delete(X[:, 0], pss.irefnode).copy()
        ## ⚠ DROP THE DUPLICATE ENDPOINT: `waveform` carries t=0 and t=T, the
        ## same point on a periodic solution, and keeping both biases every bin.
        i_probe = X[row, :-1]
        N = i_probe.shape[0]
        n = np.arange(N)
        ks = self.tones if upto is None else list(range(1, int(upto) + 1))
        I = np.array([(2.0 / N) * np.sum(i_probe * np.exp(-2j * np.pi * k * n / N))
                      for k in ks])
        return I, pss

    def probe_current(self, A, f):
        """The FUNDAMENTAL phasor of the probe current. `I1 == 0` is the
        oscillation condition.  Returns `(I1, pss)`."""
        amps = [float(A)] + [0.0] * (len(self.tones) - 1)
        phases = [self.phase] * len(self.tones)
        I, pss = self._spectrum(f, amps, phases)
        return complex(I[0]), pss

    def even_harmonic_content(self, f, amps=None, phases=None, upto=6):
        """`|I_even| / |I_odd|` on the probe current -- the ONLY basis on which
        even tones may be dropped.

        ⚠⚠ PRUNING THE EVEN HARMONICS DOES **NOT** HOLD IN GENERAL, and the
        failure is silent: dropping a tone that is really there removes both an
        unknown and the residual row constraining it, so the solve converges to
        the wrong waveform.  Van der Pol is HALF-WAVE SYMMETRIC and its even
        content is 7.1e-16; add an even term ``beta u^2`` to the same
        nonlinearity and it is not::

            beta    H2/H1       H3/H1       H4/H1
            0.00    7.080e-16   1.168e-01   2.710e-16   <- symmetric
            0.05    2.649e-02   1.162e-01   9.239e-03
            0.20    1.058e-01   1.068e-01   3.600e-02   <- H2 EQUALS H3
            0.50    2.613e-01   5.994e-02   7.610e-02   <- H2 is 4x H3

        At ``beta = 0.5`` pruning would discard the LARGEST correction after
        the fundamental, and ``beta = 0.05`` already gives 2.6% -- there is no
        margin to judge by eye.  So this measures rather than assumes, and a
        caller passing `tones=[1, 3, 5]` should check it first.
        """
        n = len(self.tones)
        if amps is None:
            amps = [1.0] + [0.0] * (n - 1)
        if phases is None:
            phases = [self.phase] * n
        I, _pss = self._spectrum(f, amps, phases, upto=upto)
        mag = np.abs(I)
        even = mag[1::2]
        odd = mag[0::2]
        ref = float(np.max(odd)) if odd.size else 0.0
        return (float(np.max(even)) / ref if ref > 0 else float('inf')), mag

    def degenerate_placement(self, A, f, tol=1e-10):
        """`(is_degenerate, periodicity_error, converged)` for this placement.

        ⚠ A placement leaving some state's DC level unconstrained gives a
        SINGULAR shooting Jacobian: every member of a one-parameter family
        satisfies periodicity, so the solve cannot converge even though what it
        returns is already a periodic solution.  The signature is exactly that
        pairing -- a tiny periodicity error with `converged = False`.
        """
        import warnings as _w
        n = len(self.tones)
        cir = self._build(f, [float(A)] + [0.0] * (n - 1),
                          [self.phase] * n)
        pss = PSS(cir, method=self.method, reltol=self.reltol)
        T = 1.0 / float(f)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / self.npts,
                      maxiterations=self.maxiterations)
        X = np.asarray(pss.waveform[1], dtype=float)
        perr = float(np.max(np.abs(X[:, -1] - X[:, 0])))
        return (perr < tol and not pss.converged), perr, bool(pss.converged)

    ## ---- solves ----------------------------------------------------------

    def _pack(self, f, amps, phases):
        """Unknowns: `f`, `A_1`, then `(A_k, phi_k)` for the rest.

        `phi_1` is NOT an unknown: it is the time origin, and carrying it would
        make the system singular along the trivial time-shift direction --
        the same marginal mode a phase condition removes in an autonomous
        solve.
        """
        z = [float(f), float(amps[0])]
        for j in range(1, len(self.tones)):
            z += [float(amps[j]), float(phases[j])]
        return np.array(z, dtype=float)

    def _unpack(self, z):
        f = float(z[0])
        amps = [float(z[1])]
        phases = [self.phase]
        for j in range(1, len(self.tones)):
            amps.append(float(z[2 * j]))
            phases.append(float(z[2 * j + 1]))
        return f, amps, phases

    def solve_multitone(self, freq0, amps0, phases0=None, tol=1e-8,
                        maxiter=20, rel_step=1e-5, use_pac=False):
        """Newton driving the probe current to zero at EVERY configured tone.

        `2K` unknowns against `2K` residuals, so the system is square.  Returns
        `(f, amps, phases, info)`.

        Measured on van der Pol (`mu = 1`, autonomous `f = 0.150229`)::

            K   f          df/f         note
            1   0.159134   +5.93e-02
            2   0.159134   +5.93e-02    A_2 = 2.3e-13 -- no change at all
            3   0.150172   -3.77e-04    157x better

        ⚠ K=2 buys NOTHING here because the even harmonics do not exist on a
        half-wave-symmetric circuit -- NOT because two tones cannot help.  See
        :meth:`even_harmonic_content` before concluding the same elsewhere.
        """
        n = len(self.tones)
        if phases0 is None:
            phases0 = [self.phase] * n
        z = self._pack(freq0, amps0, phases0)
        hist = []
        r = None
        for it in range(int(maxiter)):
            f, amps, phases = self._unpack(z)
            I, _p = self._spectrum(f, amps, phases)
            r = np.concatenate([[c.real, c.imag] for c in I])
            hist.append((f, list(amps), float(np.linalg.norm(r))))
            if np.linalg.norm(r) < tol:
                return f, amps, phases, {
                    'iterations': it, 'residual': float(np.linalg.norm(r)),
                    'history': hist, 'converged': True,
                    'evaluations': self.evaluations}
            J = np.zeros((2 * n, z.shape[0]))
            if use_pac:
                ## ⚠ The VOLTAGE block from K linear PAC solves; the FREQUENCY
                ## column stays a finite difference, because `df` moves every
                ## tone at once and is not a small-signal excitation about the
                ## same operating point.  Validated on the FIRST iteration only
                ## -- the conventions do not change as the Newton walks, and
                ## paying the check every iteration would give back the saving.
                Jv, vinfo = self.pac_jacobian(f, amps, phases,
                                              validate=(it == 0))
                ## CHAIN RULE, NOT A POSITIONAL COPY.  `pac_jacobian` returns
                ## `d(Re I, Im I)/d(Re V, Im V)`; the unknowns here are
                ## `(A_k, phi_k)` with `phi` in DEGREES.  Copying the columns
                ## across positionally feeds the Newton a Jacobian for the
                ## WRONG VARIABLES -- it diverged to f = 0.0348 against 0.1502
                ## while `pac_jacobian` itself validated at 1e-04, which is how
                ## a correct derivative and a broken solve coexisted.
                ## With `V_k = A_k exp(j psi_k)`, `psi_k = (phi_k - phase0)*pi/180`:
                ##     dV/dA   = (cos psi, sin psi)
                ##     dV/dphi = A (pi/180) (-sin psi, cos psi)
                dzf = rel_step * max(abs(z[0]), 1e-3)
                zz = z.copy()
                zz[0] += dzf
                f2, a2, p2 = self._unpack(zz)
                I2, _ = self._spectrum(f2, a2, p2)
                r2 = np.concatenate([[c.real, c.imag] for c in I2])
                J[:, 0] = (r2 - r) / dzf
                rad = np.pi / 180.0
                for jx in range(n):
                    psi = (phases[jx] - self.phase) * rad
                    cA, sA = np.cos(psi), np.sin(psi)
                    col_re = Jv[:, 2 * jx]
                    col_im = Jv[:, 2 * jx + 1]
                    dA_col = cA * col_re + sA * col_im
                    dP_col = amps[jx] * rad * (-sA * col_re + cA * col_im)
                    if jx == 0:
                        J[:, 1] = dA_col          # phi_1 is the time origin
                    else:
                        J[:, 2 * jx] = dA_col
                        J[:, 2 * jx + 1] = dP_col
            else:
                for j in range(z.shape[0]):
                    dz = rel_step * max(abs(z[j]), 1e-3)
                    zz = z.copy()
                    zz[j] += dz
                    f2, a2, p2 = self._unpack(zz)
                    I2, _ = self._spectrum(f2, a2, p2)
                    r2 = np.concatenate([[c.real, c.imag] for c in I2])
                    J[:, j] = (r2 - r) / dz
            step, _res, _rank, _sv = np.linalg.lstsq(J, -r, rcond=None)
            z = z + step
        f, amps, phases = self._unpack(z)
        return f, amps, phases, {
            'iterations': maxiter,
            'residual': float(np.linalg.norm(r)) if r is not None else None,
            'history': hist, 'converged': False,
            'evaluations': self.evaluations}

    def solve(self, amp0, freq0, tol=1e-9, maxiter=20, damp=1.0,
              rel_step=1e-4):
        """Single-tone Newton on `(A, f)`.  Returns `(A, f, info)`.

        Kept as its own entry point because the one-tone case is the cheap
        screen -- three solves an iteration -- and because its return shape is
        two scalars rather than the vectors :meth:`solve_multitone` returns.
        """
        A, f = float(amp0), float(freq0)
        hist = []
        r = np.array([np.inf, np.inf])
        for it in range(int(maxiter)):
            I0, _p = self.probe_current(A, f)
            r = np.array([I0.real, I0.imag], dtype=float)
            hist.append((A, f, float(np.linalg.norm(r))))
            if np.linalg.norm(r) < tol:
                return A, f, {'iterations': it,
                              'residual': float(np.linalg.norm(r)),
                              'history': hist, 'converged': True,
                              'evaluations': self.evaluations}
            dA = rel_step * max(abs(A), 1e-12)
            df = rel_step * max(abs(f), 1e-12)
            IA, _ = self.probe_current(A + dA, f)
            IF, _ = self.probe_current(A, f + df)
            J = np.array([[(IA.real - I0.real) / dA, (IF.real - I0.real) / df],
                          [(IA.imag - I0.imag) / dA, (IF.imag - I0.imag) / df]])
            if abs(np.linalg.det(J)) < 1e-300:
                raise np.linalg.LinAlgError(
                    'ProbeShooting: the 2x2 probe Jacobian is singular at '
                    'A=%.6g f=%.6g. Either the probe cannot see the '
                    'oscillation, or the placement leaves a state undetermined '
                    '-- check `degenerate_placement`.' % (A, f))
            step = np.linalg.solve(J, -r)
            A += damp * step[0]
            f += damp * step[1]
        return A, f, {'iterations': maxiter,
                      'residual': float(np.linalg.norm(r)),
                      'history': hist, 'converged': False,
                      'evaluations': self.evaluations}

    #: excitation offset as a fraction of `f`, so the folded sideband pair
    #: lands at distinguishable frequencies -- see `_pac_response`.
    _pac_delta = 1e-4

    def _pac_response(self, pss, cir, row, f, excite_harmonic, want_harmonics):
        """`{m: dI_m}` from ONE linear PAC solve exciting harmonic `j`.

        ⚠⚠ TWO CONVENTION FACTORS AND ONE INDEXING TRAP, all three pinned
        against a circuit whose answer is analytic (a resistor across the
        probe, where `dI/dV = 1/R` exactly) rather than against the finite
        difference this is meant to replace.  Calibrating against FD would
        make the agreement circular and would absorb a genuine sideband-index
        error into the fitted constant.

        * **The frequency list carries DUPLICATES.**  `0.159155` appears twice
          in the returned sweep, one entry near zero and one carrying the
          response; `argmin(|fs - target|)` picks whichever comes first and it
          was the wrong one, reading 5.9e-21 where the answer is 1e-3.  So the
          entry is chosen by LARGEST RESPONSE among those at the target
          frequency, not by proximity alone.
        * **A factor of two**: this file's probe spectrum uses the peak-amplitude
          convention `(2/N) sum(...)`; PAC returns a phasor.
        * **A 90 degree rotation that is OURS, not PAC's**: `VS` builds its AC
          phasor as `vac * exp(j*phase)`, and the operating-point probe sets
          `phase = 90` to make the forcing a cosine -- so the same parameter
          rotates the AC excitation.  Dividing by the excitation phasor removes
          it, which is why the excitation is read from the circuit rather than
          assumed to be 1.
        """
        tk = cir.toolkit
        pac = PAC(cir, toolkit=tk)
        ## ⚠⚠⚠ EXCITE SLIGHTLY OFF THE HARMONIC, WHICH REMOVES THE AMBIGUITY
        ## INSTEAD OF GUESSING IT.  Exciting exactly at `j*f0` sends TWO
        ## sidebands to the same absolute output frequency -- `k = m - j` and
        ## `k = -m - j` -- and `PAC.solve` returns absolute frequencies with the
        ## sideband index folded away, so the two arrive in an order that is not
        ## stable.  Ordering them by magnitude worked AT THE SOLUTION and failed
        ## away from it (validation 1.6e-05 at the solved amplitudes, 1.763 at
        ## the Newton's starting point), which is the kind of heuristic that
        ## passes a gate and then fails in use.
        ##
        ## With the excitation at `j*f0 + delta` every output lands at
        ## `(j+k)*f0 + delta`, all distinct, so the wanted term is simply the
        ## one nearest `m*f0 + delta` and no ordering rule is needed.  `delta`
        ## is small enough to be a negligible perturbation of the response and
        ## large enough to separate the pair.
        delta = self._pac_delta * float(f)
        f_ex = float(excite_harmonic) * float(f) + delta
        res = pac.solve(pss, freqs=np.array([f_ex]))
        fs = np.asarray(res.sweep_values, dtype=float)
        X = np.asarray(res.x)
        ## the excitation phasor actually applied, read from the circuit
        (u_ac,) = remove_row_col((cir.u(0, analysis='ac'),), pss.irefnode,
                                 tk)
        u_ac = np.asarray(u_ac, dtype=complex).ravel()
        scale = u_ac[np.argmax(np.abs(u_ac))]
        out = {}
        for m in want_harmonics:
            ## ⚠⚠ BOTH TERMS ARE PHYSICAL, AND THE OFFSET IS WHAT NAMES THEM.
            ## Exciting at `j*f0 + delta`, the DIRECT sideband `k = m - j`
            ## lands at `m*f0 + delta`, and the IMAGE `k = -m - j` lands at
            ## `-m*f0 + delta`, which PAC folds onto `m*f0 - delta` and
            ## conjugates.  Taking only the direct one halves the answer
            ## (measured: a uniform 0.5 at every K); taking both without being
            ## able to tell them apart is what forced the magnitude-ordering
            ## heuristic that passed at the solution and failed at the Newton's
            ## start.  With the offset they are identified by FREQUENCY, so the
            ## rule is derived rather than guessed.
            base = float(m) * float(f)
            i_dir = np.where(np.abs(fs - (base + delta)) < 0.25 * delta)[0]
            i_img = np.where(np.abs(fs - (base - delta)) < 0.25 * delta)[0]
            if i_dir.size == 0:
                continue
            d_val = complex(X[row, i_dir[int(np.argmin(np.abs(
                fs[i_dir] - (base + delta))))]])
            g_val = (complex(X[row, i_img[int(np.argmin(np.abs(
                fs[i_img] - (base - delta))))]]) if i_img.size else 0.0 + 0.0j)
            out[m] = -(d_val - g_val) / scale
        return out

    def pac_jacobian(self, f, amps, phases, validate=True, rel_step=1e-5,
                     rtol=5e-2):
        """`dI/dV` from K LINEAR PAC solves instead of 2K nonlinear ones.

        The finite-difference Jacobian recomputes, by `2K` full PSS solves, a
        quantity that is a PERIODIC SMALL-SIGNAL response about the operating
        point the base solve already produced.  One PSS plus K PAC solves gets
        the voltage block, which is where the cost is.

        Returns `(J, info)` with `J` the real `2K x 2K` block
        `d(Re I, Im I)/d(Re V, Im V)`.

        ⚠⚠⚠ **NOT SHIPPED-READY: THE NORMALISATION IS INCOMPLETE, AND THE
        DEFAULT VALIDATION CORRECTLY REFUSES.**  The route is confirmed viable
        -- PAC returns the right quantity, verified against a resistor where
        `dI/dV = 1/R` analytically -- and two of the three discrepancies are
        pinned and removed.  A third is not:

            fixture                 PAC / FD after normalisation
            resistor (linear)       -1.0        (sign only)
            van der Pol (K=1)       -0.499999   (sign AND a factor 2)

        **A constant that differs between two circuits is not a convention,
        it is a missing term**, so the remaining factor is NOT applied by
        fitting it -- that would make the "independent" Jacobian a fit to the
        finite difference it replaces, hide any sideband-index error inside the
        fitted constant, and reproduce exactly the circular-verification
        failure this file already records for the PPV.

        Until it is derived, `validate=True` raises on real circuits and the
        finite-difference Jacobian in :meth:`solve_multitone` remains the
        shipped path.  What IS established and reusable:

        * the response is present and correct (1/R recovered exactly);
        * the 90 degree rotation is OURS -- `VS` builds its AC phasor as
          `vac * exp(j*phase)` and the operating-point probe sets `phase = 90`
          -- and dividing by the excitation phasor removes it, measured: the
          residual ratio is real, not imaginary;
        * PAC's frequency sweep contains DUPLICATE entries at the same
          frequency, one near zero and one carrying the response, so
          `argmin(|fs - target|)` reads 5.9e-21 where the answer is 1e-3.

        ⚠⚠ `validate=True` CHECKS ONE COLUMN AGAINST THE FINITE DIFFERENCE AND
        RAISES ON DISAGREEMENT, and it is on by default deliberately.  A wrong
        Jacobian does not announce itself: the Newton still converges, to the
        wrong orbit -- the same silent failure that even-harmonic pruning
        produces, which this file already has a falsifier for.  The check costs
        ONE extra solve, amortised over the whole solve, and it exercises the
        conventions in `_pac_response` on the circuit actually in hand rather
        than on the resistor they were derived from.

        ⚠ The OFF-DIAGONAL entries `dI_m/dV_j` with `m != j` are the ones that
        exercise the sideband map `k = m - j` and PAC's negative-frequency
        conjugation.  With `harmonics=1` there are none, so a passing K=1
        validation says NOTHING about the index map -- validate at K >= 2
        before trusting a multitone Jacobian.
        """
        import warnings as _w
        n = len(self.tones)
        ## ⚠⚠ ONE PSS SOLVE FOR EVERY COLUMN.  `vac` is read ONLY under
        ## `analysis='ac'` -- it does not enter the transient residual, so it
        ## cannot move the periodic operating point.  Building a fresh circuit
        ## and re-solving the PSS per column therefore recomputed the SAME
        ## orbit K times, which is the whole cost this method exists to avoid:
        ## it made the "cheap" Jacobian K nonlinear solves plus K linear ones,
        ## against the finite difference's 2K.  Solve once, then walk `vac`
        ## across the probes and take K LINEAR PAC solves against that one
        ## operating point.
        cir = self._build(f, amps, phases)
        row = self._probe_row(cir)
        pss = PSS(cir, method=self.method, reltol=self.reltol)
        T = 1.0 / float(f)
        with _w.catch_warnings():
            _w.simplefilter('ignore')
            pss.solve(period=T, timestep=T / self.npts,
                      maxiterations=self.maxiterations)
        if not pss.converged:
            raise RuntimeError(
                'pac_jacobian: the periodic operating point did not converge, '
                'so there is nothing to linearise about.')
        J = np.zeros((2 * n, 2 * n))
        cols = {}
        for jx, k in enumerate(self.tones):
            ## excite exactly one probe; the others must be silent, and `vac`
            ## DEFAULTS TO 1 on a VS, which is why every one is set explicitly
            for jj in range(n):
                cir['__probe%d' % jj].ipar.vac = 1.0 if jj == jx else 0.0
            cir.update_iparv()
            resp = self._pac_response(pss, cir, row, f, k, self.tones)
            cols[jx] = resp
            for ix, m in enumerate(self.tones):
                d = resp.get(m, 0.0 + 0.0j)
                J[2 * ix, 2 * jx] = d.real
                J[2 * ix + 1, 2 * jx] = d.imag
                J[2 * ix, 2 * jx + 1] = -d.imag
                J[2 * ix + 1, 2 * jx + 1] = d.real
        info = {'columns': cols, 'validated': False}
        if validate:
            I0, _ = self._spectrum(f, amps, phases), None
            I0 = I0[0]
            dA = rel_step * max(abs(amps[0]), 1e-12)
            a2 = list(amps)
            a2[0] += dA
            I1, _p = self._spectrum(f, a2, phases)
            fd = (I1 - I0) / dA
            pac_col = np.array([J[2 * ix, 0] + 1j * J[2 * ix + 1, 0]
                                for ix in range(n)])
            num = np.linalg.norm(pac_col - fd)
            den = max(np.linalg.norm(fd), 1e-300)
            info['validation_reldiff'] = float(num / den)
            info['validated'] = bool(num / den < rtol)
            if not info['validated']:
                raise ValueError(
                    'pac_jacobian: the PAC column disagrees with the finite '
                    'difference by %.3e (relative), above rtol=%.2e. PAC gave '
                    '%s against FD %s. Do NOT use this Jacobian -- a wrong one '
                    'converges silently to the wrong orbit.'
                    % (num / den, rtol,
                       np.array2string(pac_col, precision=5),
                       np.array2string(fd, precision=5)))
        return J, info

    def power_flow(self, A, f, rel_step=1e-4):
        """The 2x2 power-flow screen `P = dv_R dy_R + dv_I dy_I`.

        ⚠⚠ ONE-DIRECTIONAL, AND THAT IS THE AUTHORS' OWN STATEMENT.  Only
        `P > 0 => unstable` is proven; they say the converse and
        `P < 0 => stable` "have not been proven", only tested.  So this is a
        cheap INSTABILITY DETECTOR and never a replacement for
        `_spectral_report`: a non-positive `P` means NOT DETECTED, not stable.

        Returns `(P_max, info)`, `P_max` being the largest value over
        perturbation directions -- the largest eigenvalue of the symmetric part
        of `dY/dV`.  It skips the system Jacobian's eigenvalues, which is the
        whole point of the construction.
        """
        I0, _ = self.probe_current(A, f)
        Y0 = I0 / A
        dv = rel_step * max(abs(A), 1e-12)
        Ir, _ = self.probe_current(A + dv, f)
        Yr = Ir / (A + dv)
        ph = dv / max(abs(A), 1e-12)
        Ii, _ = self.probe_current(A, f * (1.0 + ph / (2.0 * np.pi)))
        Yi = Ii / A
        dYdV = np.array([[(Yr.real - Y0.real) / dv, (Yi.real - Y0.real) / dv],
                         [(Yr.imag - Y0.imag) / dv, (Yi.imag - Y0.imag) / dv]])
        sym = 0.5 * (dYdV + dYdV.T)
        eig = np.linalg.eigvalsh(sym)
        return float(np.max(eig)), {'dYdV': dYdV, 'symmetric_part': sym,
                                    'eigenvalues': eig, 'Y0': Y0,
                                    'unstable': bool(np.max(eig) > 0.0)}
