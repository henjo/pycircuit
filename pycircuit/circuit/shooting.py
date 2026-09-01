from pycircuit.post import InternalResultDict
from .circuit import gnd
from pycircuit.circuit.analysis import *
from copy import copy
import warnings
import pycircuit.circuit.analysis as analysis
import numpy as np

def freq_analysis(x, t, rms = True, axis=-1, freqoffset = 0):
    """Return dft of equidistant sampled signal x"""
    
    npoints = np.size(x, axis)

    dt = t[1] - t[0]

    if x.dtype in (np.cdouble, np.cdouble):
        X = np.fft.fftshift(np.fft.fft(x, axis=axis),axes=(axis,)) / npoints
        freqs = np.fft.fftshift(np.fft.fftfreq(npoints, d=dt))
    else:
        freqs = np.fft.fftfreq(npoints, d=dt)[:int(np.ceil(npoints / 2.))]
        slices = [slice(None)] * x.ndim
        slices[axis] = slice(0, len(freqs))
        X = np.fft.fft(x, axis=axis)[tuple(slices)] / npoints
        ## Fold energy from negative frequencies
        X[:,1:] *= np.sqrt(2)

    if not rms:
        X *= np.sqrt(2)

    return freqs, X

## ⚠ THE SPECTRAL RADIUS CANNOT DECIDE THIS, and trying it first is the
## instructive part.  An autonomous orbit gives an eigenvalue at exactly 1
## -- but only AT its own period, and a run at any other period reads well
## below (measured 0.9615 on the quadrature phase element at the nominal
## period against 1.000226 at the corrected one).  Worse, a merely
## lightly-damped DRIVEN circuit sits near 1 too: a Q=1000 resonator has
## `exp(-pi/Q) = 0.99686`.  So no threshold separates the two -- one
## setting misses the autonomous case where users will actually run it, the
## other fires on every high-Q filter.
##
## The distinction is structural, not spectral, and it is exact: a circuit
## is autonomous when nothing in it depends on `t`.  `u(t)` is sampled
## across the period and compared; a phase accumulator driven by a DC
## source is autonomous however energetically it oscillates, which is
## precisely the case arc 5 asks about.
AUTONOMOUS_U_TOL = 1e-12


class PSS(Analysis):
    """Periodic Steady-State using shooting Newton iterations

    The algorithm is described in [1] p65.

     1. Kenneth S. Kundert, Jacob K. White, Alberto Sangiovanni-Vincentelli
        (1990)
        Steady-State Methods for Simulating Analog and Microwave Circuits
        Kluwer Academic Publishers
        ISBN 0792390695

    **THREE CONVERGENCE CHECKS NEST HERE, AND THEY ARE NOT INTERCHANGEABLE.**

    1. the per-timestep Newton, which solves the discretised circuit
       equations at one time point;
    2. the local truncation error, which decides how far the DISCRETE
       trajectory is from the true one;
    3. the shooting Newton, which finds the periodic point of the discrete
       map.

    Two rules order them.  **(3) cannot be tighter than (1)**: the period
    map is only KNOWN to the accuracy of the per-timestep solves, so the
    shooting residual has a floor there whatever its Jacobian is.  And
    **LTE must not run per shooting iteration**: an
    adaptive grid makes the step sequence a function of `x0`, so the period
    map stops being smooth and (3) loses its quadratic rate.  Choose the
    grid once, freeze it, shoot on it.

    **(3) is a true Newton for both methods.**  It was neither, for years:
    the monodromy accumulation was missing a factor of `1/h`, and since `C`
    is singular the product collapsed to exactly zero -- so the Jacobian was
    the identity and this was successive substitution, `x0 <- phi(x0)`, on
    every circuit and without saying so.  It was found by its RATE: the
    residual fell 0.855 per iteration on a Q=20 resonator, and `exp(-pi/Q)
    = 0.8546` is that circuit's own per-period decay.  Non-convergence is
    reported now; it used to be discarded with `full_output=False`.

    Trapezoidal needed more than the factor, and Gear-2 more again.  Every
    method here writes its companion as

        iq_n = sum_k a_k q_{n-k}  +  b iq_{n-1}

    so ONE recursion differentiates all of them:

        S    = sum_{k>=1} a_k C_{n-k} Px_{n-k}  +  b Pq
        Px_n = -Jf_n^-1 S
        Pq_n = a_0 C_n Px_n + S

    Euler is `b = 0` reaching back one step, trapezoidal `b = -1` reaching
    back one, Gear-2 `b = 0` reaching back two.  The coefficients come from
    `Integrator.companion_coefficients` -- from the integrator that ACTUALLY
    ran, so an order-dropped opening step contributes its own -- rather than
    being transcribed here, which this tree has paid for three times.

    An x-only monodromy is not merely less accurate for the methods with
    memory: measured, applying the Euler form to trapezoidal converged
    SLOWER than no Jacobian at all (0.90 against 0.855).  A wrong Jacobian
    is worse than none.

    ⚠ **THE TWO FAILURE MODES ARE ORTHOGONAL, and level 3 cannot see level
    2.**  On the Q=20 resonator against a 20 V analytic peak, all three now
    converge -- and they do not agree:

        euler   5 iterations,  8.815 V   (56% low)
        gear2   7 iterations, 19.766 V   (1.2% low)
        trap    6 iterations, 19.990 V   (0.05% low)

    That is each method's own numerical damping, invisible to levels 1 and
    3, and it is why the LTE report in the recorded scope below is worth
    more than it looks: a converged shooting solve is not by itself
    evidence of a correct answer.

    (2) DOES NOT EXIST HERE, AND UNDER A SHOOTING METHOD IT CANNOT BE A
    CONTROLLER.  The grid is a fixed `linspace`; under `fixed_timestep` both
    transient backends skip the LTE verdict outright, keeping only the
    order drop that protects the integrator across a breakpoint.  Nor is
    that a wiring gap to be closed: if the step sequence adapts to `x0`
    then phi is a DIFFERENT discrete map for each `x0`, so it is not smooth
    in `x0` and the accumulated monodromy is the derivative of a
    neighbouring problem.  Freezing the grid is what makes (3) a Newton.

    So (2) changes kind here.  The estimator is still computable on a
    frozen grid, and what it measures is not convergence but ACCURACY:
    which of the three levels is limiting the answer.  That question is
    live rather than academic -- on a Q=20 resonator with `method='euler'`
    the shooting solve converges completely (5 iterations, residual
    3.9e-05) and lands at 8.815 V against a 20 V analytic peak.  The answer
    is 56% low for a reason (1) and (3) cannot see, and nothing currently
    says so.

    RECORDED SCOPE, in order, none of it planned work yet:

      4. LTE as a REPORT -- evaluate the estimator on the converged
         periodic solution and say whether the run is discretisation-
         limited.  Needs only the estimator exposed; catches the case above.
      5. LTE-CHOSEN grid -- pick the step sequence from an adaptive run and
         freeze it, refining BETWEEN shooting solves.  The grid still never
         moves inside one, so (3) stays exact.  Blocked on `Transient`
         accepting a non-uniform grid; `fixed_timestep` is uniform-only.
         Note this is the better structure anyway: a transient adapts
         because it cannot see the future, while PSS re-solves the same
         interval repeatedly and can therefore choose its grid once, well.
      6. MATRIX-FREE variational shooting (Telichevesky, Kundert & White,
         DAC 1995) -- the only structure that permits per-iteration
         adaptivity, because `M v` is obtained by integrating the
         variational system along the stored trajectory on that iteration's
         own grid.  The outer solve becomes an INEXACT Newton, phi shifting
         slightly between iterations.  A rewrite, not an increment on the
         above, and it should not leak into one.

    Driving `Transient` at `fixed_timestep=True` -- the next planned step --
    buys one integrator definition, the limiting/PCNR machinery, breakpoints
    and the order drop.  It does NOT buy (2); saying otherwise was this
    docstring's own earlier error.

    `reltol`, `iabstol` and `vabstol` mean exactly what they mean on
    `Transient` -- the tolerances of the TRANSIENT solution, applied to the
    per-timestep Newton, with the two absolute floors applied PER UNKNOWN in
    both flavours by `analysis.newton_tolerance_vectors`, the single
    definition all three analyses read.  Nothing here rescales them.

    The shooting criterion is expressed against that one: `steadyratio`
    (>= 1, default 1) multiplies it, so by default the shooting solve is
    held to the SAME relative tolerance as the transient, and raising it
    buys fewer shooting iterations for a looser periodic steady state.
    """

    parameters = Analysis.parameters + \
        [Parameter(name='analysis', desc='Analysis name',
                   ## Sources supply their time-domain waveform only for an
                   ## analysis name in timedomain_analyses (('dc','tran')); the
                   ## old default 'PSS' matched nothing, so cir.u(t) returned 0
                   ## and the whole shooting solve had no excitation.
                   default='tran'),
         Parameter(name='reltol', 
                   desc='Relative tolerance', unit='', 
                   default=1e-4),
         Parameter(name='iabstol', 
                   desc='Absolute current error tolerance', unit='A', 
                   default=1e-12),
         Parameter(name='vabstol', 
                   desc='Absolute voltage error tolerance', unit='V', 
                   default=1e-12),
         Parameter(name='maxiter', 
                   desc='Maximum number of iterations', unit='', 
                   default=100),
         ## `reltol` MEANS THE SAME THING IN EVERY ANALYSIS: the relative
         ## tolerance of the transient solution.  It is applied to the
         ## per-timestep Newton here exactly as `Transient` applies it, and
         ## nothing rescales it.
         ##
         ## `steadyratio` is how the SHOOTING criterion is expressed relative
         ## to it: shooting reltol = reltol * steadyratio, with 1 meaning the
         ## two are equal.  It is >= 1 because the period map is only KNOWN
         ## to the accuracy of the inner solves, so asking the shooting
         ## residual to beat that is asking it to resolve its own noise --
         ## refused rather than silently accepted.  Raise it to accept a
         ## looser periodic steady state for fewer shooting iterations.
         Parameter(name='steadyratio',
                   desc='Shooting tolerance as a multiple of reltol (>= 1); '
                        '1 holds the shooting solve to the same relative '
                        'tolerance as the transient, larger relaxes it',
                   unit='', default=1.0),
         Parameter(name='method',
                   desc="Integration method for the inner transient: 'trap' "
                        "(default), 'euler', or 'gear' (BDF-2)",
                   unit='',
                   default="trap")]        

    
    def __init__(self, cir, toolkit=None, irefnode=None, **kvargs):
        self.parameters = super(PSS, self).parameters + self.parameters
        super(PSS, self).__init__(cir, **kvargs)
        ## The reference row is fixed for the analysis, and both the shooting
        ## loop and the Transient this drives need it.  It was recomputed in
        ## every method from a `refnode` argument that no caller ever varied.
        self.irefnode = self.cir.get_node_index(
            gnd if irefnode is None else irefnode)
        self._tran = None

    def _is_autonomous(self, times):
        """True when nothing in the circuit depends on `t`.

        Exact where a spectral test is not: see `AUTONOMOUS_U_TOL`.  The
        source vector is sampled across the period and compared with the
        first sample; a circuit driven only by DC -- a VCO macromodel, a
        phase accumulator, an LC or ring oscillator -- has a constant `u`
        and a one-parameter family of periodic solutions, which is what
        makes fixed-period shooting ill-posed for it.
        """
        u0 = np.asarray(self.cir.u(times[0], analysis=self.par.analysis),
                        dtype=float)
        scale = max(float(np.max(np.abs(u0))), 1.0)
        for t in times[1::max(1, len(times) // 8)]:
            u = np.asarray(self.cir.u(t, analysis=self.par.analysis),
                           dtype=float)
            if float(np.max(np.abs(u - u0))) > AUTONOMOUS_U_TOL * scale:
                return False
        return True

    def _transient(self):
        """The `Transient` this analysis integrates with.

        PSS used to carry its OWN transcription of one integrator step --
        the third in the tree, after `Transient` and `JAXTransient` -- and
        it had already cost the two defects its docstring records: `method`
        declared and never read, and a companion current fed back from the
        iterate before the converged one.  Driving the real thing removes
        the copy and brings what came with it: the limiting machinery, PCNR,
        breakpoint order drops, and the continuation rescue.

        The tolerances are handed over unchanged, which is the point of
        `newton_tolerance_vectors`: `reltol`/`iabstol`/`vabstol` mean the
        same thing on both sides, so passing them through is a no-op in
        meaning.
        """
        if getattr(self, '_tran', None) is None:
            from pycircuit.circuit.transient import Transient
            from pycircuit.circuit.integrator import (EulerIntegrator,
                                                      TrapezoidalIntegrator,
                                                      Gear2Integrator)
            ## ⚠ A MAPPING, not an if/else on 'euler'.  Written as
            ## `EulerIntegrator() if method == 'euler' else Trapezoidal...`
            ## it silently ran trapezoidal for every other name -- caught
            ## while adding 'gear', which produced numbers identical to
            ## trap's to the last digit.  This class has already paid once
            ## for a `method` that selected nothing; a dict raises KeyError
            ## on a name nobody wired.
            integ = {'euler': EulerIntegrator,
                     'trap': TrapezoidalIntegrator,
                     'trapezoidal': TrapezoidalIntegrator,
                     'gear': Gear2Integrator,
                     'gear2': Gear2Integrator}[self.par.method]()
            self._tran = Transient(
                self.cir, toolkit=self.toolkit, integrator=integ,
                reltol=self.par.reltol, iabstol=self.par.iabstol,
                vabstol=self.par.vabstol, maxiter=self.par.maxiter,
                analysis=self.par.analysis)
            self._tran.irefnode = self.irefnode
        return self._tran

    def _begin_period(self, x_reduced):
        """Start one traversal of the period from a clean integrator state.

        Every shooting iteration re-integrates the SAME interval from its own
        `x0`, so "begin a run" happens once per iteration here, not once per
        analysis.  Without the reset, iteration k+1 would inherit the ring
        buffers iteration k ended with and the period map would depend on
        which iteration it was -- phi must be a function of `x0` alone or the
        monodromy is the derivative of something else.
        """
        tr = self._transient()
        tr._begin_run(self._insert_refnode(x_reduced), self.cir.n)
        return tr

    def _insert_refnode(self, x):
        return self.toolkit.concatenate(
            (x[:self.irefnode], self.toolkit.array([0.0]), x[self.irefnode:]))

    def solve_timestep(self, x0, t, dt, refnode=gnd, iq_last=None):
        """One timestep of the inner transient, taken by `Transient`.

        This used to be a private transcription of one integrator step --
        the third in the tree -- and it had already cost two defects that
        its own comments recorded: `method` was declared and never read, so
        PSS was backward-Euler only, and the companion current fed back to
        the next step belonged to the iterate BEFORE the converged one.
        Both are structurally impossible now: the integrator is an
        `Integrator` object driven by `Transient.get_diff`, and the
        companion current is the one that class stores at its own converged
        point.

        What came with the change, none of which the copy had: the limiting
        machinery (measured -- a rectifier whose diode never turned on, so
        the non-conducting solution was returned as a converged periodic
        steady state), PCNR when the circuit and Parameters ask for it,
        breakpoint order drops, and the continuation rescue.

        `dt` is imposed by the caller: PSS owns the grid, which is what
        keeps the period map a smooth function of `x0`.  `iq_last` is
        retained in the signature for callers that pass it, but the
        companion history now lives in the `Transient`'s own ring buffers,
        rolled here through `_push_history`.

        The measured cost of backward Euler on a limit cycle is unchanged
        and still the reason `method` matters -- it damps exactly what PSS
        exists to find:

            steps/period    PSS peak    fraction of analytic
                      20      2.63 V       13.2%
                      50      5.61 V       28.1%
                     100      8.81 V       44.1%
                     200     12.20 V       61.0%
        """
        toolkit = self.toolkit
        irefnode = self.irefnode
        tr = self._transient()

        ## ONE INTEGRATOR STEP, taken by the class that owns the definition.
        ## `Transient.solve_timestep` applies the chosen integrator through
        ## `get_diff` (so `method` selects something because the integrator
        ## object does), the limiting machinery, PCNR when asked for, and the
        ## continuation rescue.  None of that existed on the copy this
        ## replaced.
        tr._dt = dt
        x_full, J_full = self._transient_step(tr, x0, t)

        ## The history advance is the accept path's, called rather than
        ## copied -- and `_dt_last` must roll AFTER the step, because
        ## `get_diff` read it as `h_last` while solving.
        tr._push_history(x_full)
        tr._dt_last2 = tr._dt_last
        tr._dt_last = dt
        tr._is_first_step = False
        tr._no_history = False

        ## Reduced-system views for the shooting Jacobian.  `_Geq` is the
        ## companion conductance the step actually used, which is the factor
        ## the monodromy needs; `_iq` is kept for the caller's own bookkeeping
        ## as before.
        (self._Jf, self._Geq, self._C) = remove_row_col(
            (J_full, tr._Geq, tr._Cmat), irefnode, toolkit)
        ## The coefficients of the integrator that ACTUALLY ran this step --
        ## an order drop on the opening step reports Euler's, which is what
        ## the propagation must use for that step.
        self._coeffs = tr._companion_coeffs
        self._iq = tr._iq

        x = toolkit.concatenate((x_full[:irefnode], x_full[irefnode + 1:]))
        return x

    def _transient_step(self, tr, x0_reduced, t):
        """`Transient.solve_timestep` on the FULL vector, returning
        ``(x, J)``.  PSS works on the reduced system throughout; this is the
        one place the two conventions meet."""
        x_full = self._insert_refnode(x0_reduced)
        x, _feval, J, _f = tr.solve_timestep(x_full, t)
        return x, J


    def solve(self, refnode=gnd, period=1e-3, x0=None, timestep=1e-6, 
              maxiterations=20):
        self.period = period
        toolkit = self.toolkit

        irefnode=self.cir.get_node_index(refnode)
        n = self.cir.n
        dt = timestep
        if x0 is None:
            x = toolkit.zeros(n-1) #currently without reference node !
        else:
            x = x0 # reference node not included !

        #create vector with timepoints and a more fitting dt
        times,dt=toolkit.linspace(0,period,num=int(period/dt),endpoint=True,
                             retstep=True)
        alpha = 1

        ## Resolved here as well as in `solve_timestep`, because the SHOOTING
        ## Jacobian depends on which integrator the inner steps used.
        method = getattr(self.par, 'method', 'euler')
        if method not in ('euler', 'trap', 'trapezoidal', 'gear', 'gear2'):
            raise ValueError(
                "method must be 'euler', 'trap' or 'gear', not %r" % (method,))

        ## THE SHOOTING JACOBIAN FOLLOWS THE INTEGRATOR'S OWN COEFFICIENTS.
        ##
        ## Backward Euler's per-step sensitivity is
        ##     dx_n/dx_{n-1} = Jf_n^-1 * C(x_{n-1})/h
        ## -- the COMPANION CONDUCTANCE at the previous point, not the raw
        ## capacitance matrix.  The `/h` was missing, and because C is
        ## singular the accumulated product collapsed to EXACTLY ZERO: the
        ## Jacobian handed to fsolve was `I - 0 = I`, so the "shooting
        ## Newton" was plain successive substitution `x0 <- phi(x0)`.  That
        ## converges at the circuit's own per-period decay -- measured on a
        ## Q=20 resonator as 0.855 per iteration against exp(-pi/Q) = 0.8546,
        ## which is how it was found -- and it never reached fsolve's
        ## tolerance, on any circuit, silently.  With the companion
        ## conductance the same resonator converges in FIVE iterations
        ## (2.64 -> 2.6e-2 -> 1.0e-2 -> 1.9e-4 -> 3.9e-5).
        ##
        ## TRAPEZOIDAL NEEDS THE AUGMENTED STATE, and now has it.  Its
        ## recursion carries `iq` as well as `x`, so the period map is a
        ## function of (x, iq); an x-only monodromy is structurally
        ## incomplete, and measured, using the Euler form for trap converged
        ## SLOWER than no Jacobian at all (0.90 against 0.855 per iteration).
        ##
        ## Differentiating the two recursions together,
        ##
        ##     iq_n = 2(q_n - q_{n-1})/h - iq_{n-1}
        ##     0    = i(x_n) + iq_n + u(t_n)
        ##
        ## gives a propagation of `d(x,iq)/dx0` that costs one extra matrix
        ## product over the Euler form:
        ##
        ##     rhs   = Geq_{n-1} Px + Pq
        ##     Px_n  = Jf_n^-1 rhs
        ##     Pq_n  = Geq_n Px_n - rhs
        ##
        ## Euler is the SAME recursion with `Pq == 0`: its companion carries
        ## no `iq_{n-1}` term, so dF/diq_{n-1} vanishes and the second row
        ## never enters.  One formula, two methods, which is why this is not
        ## a second code path.
        newton = True

        def func(x):
            ## EVERY SHOOTING ITERATION IS ITS OWN RUN.  phi must be a
            ## function of `x0` alone; if iteration k+1 inherited the ring
            ## buffers iteration k ended with, the period map would depend on
            ## which iteration it was and the monodromy would be the
            ## derivative of something else.  `_begin_run` also makes the
            ## first step of each traversal `is_first_step`, so trapezoidal
            ## falls back to Euler there -- the same restart the old
            ## `iq_last=None` produced, now for the integrator's own reason.
            self._begin_period(x)
            x = self.solve_timestep(x, times[0], dt)
            iq_last = self._iq
            x0 = copy(x)
            ## `Px` is d(x_n)/d(x0) and `Pq` is d(iq_n)/d(x0).  The first
            ## step of the period is taken by Euler (the integrator drops
            ## order on `is_first_step`), so the companion entering the loop
            ## is `(q_1 - q_0)/h` and its sensitivity to `x_1` is that step's
            ## own `Geq`.  Under Euler the row is identically zero and stays
            ## so.
            ## `Px_k` is d(x_{n-k})/d(x0) and `Pq` is d(iq_n)/d(x0); `Cs_k`
            ## is the capacitance matrix of step n-k.  Two of each, which is
            ## as far back as any method here reaches.  Both rings open
            ## seeded with the entering step, mirroring how the transient
            ## seeds `_qlast` with `q0` repeated -- at the start of a period
            ## there is no earlier point to differentiate against.
            eye = np.asarray(toolkit.eye(n - 1))
            Px = [eye, eye]
            Cs = [copy(np.asarray(self._C)), copy(np.asarray(self._C))]
            a_first, b_first = self._coeffs
            Pq = (a_first[0] * Cs[0] if b_first else
                  np.zeros((n - 1, n - 1)))

            ## Save C and transient jacobian for PAC analysis
            self.Cvec = [copy(self._C)]
            self.Jtvec = [copy(self._Jf)]
            self.times = times
            for t in times[1:]:
                x = copy(self.solve_timestep(x, t, dt, iq_last=iq_last))
                iq_last = self._iq
                self.Cvec.append(copy(self._C))
                self.Jtvec.append(copy(self._Jf))
                ## STAGE 11 -- A SOLVE, NOT AN INVERSE.
                ##
                ## This read `inv(Jf) @ C @ Jshoot`, which forms an explicit dense
                ## inverse at every timestep of every shooting iteration: at
                ## N=137, M=1000 with 20 shooting iterations that is 20,000 dense
                ## inversions plus 40,000 dense matmuls.  The quantity wanted is
                ## the solution of `Jf X = C @ Jshoot`, and asking for it directly
                ## is both faster and better conditioned -- forming an inverse
                ## squares the condition number you then multiply through.
                ##
                ## The standard result in the field (Telichevesky, Kundert & White,
                ## DAC 1995) is stronger still: matrix-free shooting needs only
                ## products with the monodromy matrix, never the matrix itself.
                ## That is a rewrite; this is the same computation, correctly
                ## expressed, and it is bit-comparable rather than merely close.
                ## ONE RECURSION FOR EVERY METHOD.  Each writes its
                ## companion as `iq_n = sum_k a_k q_{n-k} + b iq_{n-1}`, so
                ## differentiating the step gives
                ##
                ##     S    = sum_{k>=1} a_k C_{n-k} Px_{n-k} + b Pq
                ##     Px_n = -Jf_n^-1 S
                ##     Pq_n = a_0 C_n Px_n + S
                ##
                ## Euler is `b = 0` with one past term, trapezoidal `b = -1`
                ## with one, Gear-2 `b = 0` with two.  The coefficients come
                ## from the integrator that ran, so an order-dropped step
                ## contributes its own.
                alphas, b = self._coeffs
                S = b * Pq if b else np.zeros_like(Px[0])
                for k in range(1, len(alphas)):
                    S = S + alphas[k] * (Cs[k - 1] @ Px[k - 1])
                Px_new = -self.toolkit.linearsolver(np.asarray(self._Jf), S)
                C_new = np.asarray(self._C)
                Pq = alphas[0] * (C_new @ Px_new) + S
                Px = [Px_new, Px[0]]
                Cs = [copy(C_new), Cs[0]]

            residual = x0 - x

            D = np.asarray(toolkit.eye(n-1))
            ## Kept for the autonomous check after the solve: its spectrum is
            ## the only place a free period announces itself.
            self._monodromy = Px[0]
            return residual, D - alpha * Px[0]
        
        ## THE SHOOTING RESIDUAL IS IN SOLUTION UNITS, NOT KCL UNITS.
        ## `x0 - phi(x0)` is a difference of SOLUTIONS -- volts on node rows,
        ## amps on branch rows -- so its absolute floor is the `xtol` flavour
        ## (vabstol on nodes, iabstol on branches), not the residual flavour
        ## the transient's Newton uses for `i(x)`.  Getting that backwards is
        ## F6(a)'s defect, and it is easy to walk into here because the
        ## quantity is called a residual.
        _tol = analysis.newton_tolerance_vectors(
            len(self.cir.nodes), len(self.cir.branches),
            self.par.iabstol, self.par.vabstol, self.toolkit)[1]
        (_tol,) = remove_row_col((_tol,), irefnode, self.toolkit)

        ## The shooting criterion, expressed against the transient one.
        _ratio = float(self.par.steadyratio)
        if _ratio < 1.0:
            raise ValueError(
                'steadyratio must be >= 1 (got %g): the period map is only '
                'known to the accuracy of the per-timestep solves, so a '
                'shooting tolerance tighter than reltol asks the outer '
                'residual to resolve its own noise.' % _ratio)
        _shoot_reltol = self.par.reltol * _ratio
        _tol = _tol * _ratio

        ## Find periodic steady state x-vector
        x0_ss, _info, _ier, _mesg = analysis.fsolve(
            func, x, maxiter=maxiterations, reltol=_shoot_reltol,
            abstol=_tol, xtol=_tol, toolkit=self.toolkit, full_output=True)
        self.converged = (_ier == 1)
        self.shooting_iterations = maxiterations if not self.converged else None
        ## ⚠ AN AUTONOMOUS OSCILLATOR CANNOT BE SOLVED AT A FIXED PERIOD,
        ## and this is the only place it says so.
        ##
        ## A circuit whose oscillation is self-sustaining -- a VCO
        ## macromodel, an LC or ring oscillator, any phase accumulator
        ## driven by a DC source -- has a one-parameter family of periodic
        ## solutions, because rotating the starting point along the orbit
        ## gives another one.  Its monodromy therefore has an eigenvalue at
        ## exactly 1 and `I - M` is singular AT the true period.  Away from
        ## it the orbit does not close at all: measured on the quadrature
        ## phase element, the discretisation precesses by 2.1e-3 rad per
        ## cycle at 100 steps/period (falling as h^2), so the period map is
        ## a rotation by slightly less than 2*pi whose ONLY fixed point is
        ## the origin -- which is what an unseeded run returns, silently.
        ##
        ## Measured either side on that element: |eig(M)| = 0.968 and
        ## sigma_min(I-M) = 2.3e-02 at the nominal period, against
        ## |eig(M)| = 1.000226 and 1.6e-04 at the corrected one.  So there is
        ## no period at which this analysis both has a solution and an
        ## invertible Jacobian, and the answer is not a better seed: it is
        ## the autonomous formulation, which solves for the period jointly
        ## with a phase condition.  Not implemented -- but a run that
        ## returns the origin, or refuses to converge, deserves to be told
        ## why rather than left to look like a tolerance problem.
        rho = None
        if getattr(self, '_monodromy', None) is not None:
            try:
                rho = float(np.max(np.abs(np.linalg.eigvals(
                    np.asarray(self._monodromy)))))
            except np.linalg.LinAlgError:                 # pragma: no cover
                rho = None
        self.spectral_radius = rho
        self.autonomous = self._is_autonomous(times)
        if self.autonomous:
            warnings.warn(
                'PSS: no source in this circuit depends on time, so it is '
                'AUTONOMOUS (self-oscillating) and a FIXED-period shooting '
                'analysis cannot solve it. At the orbit period the starting '
                'phase is free, so the Jacobian is singular; at any other '
                'period the orbit does not close and the only periodic '
                'solution is the trivial one -- an unseeded run returns it. '
                'Solving for the period jointly with a phase condition is a '
                'different analysis and is not implemented. The monodromy '
                'spectral radius here is %s; the result below is not a '
                'periodic steady state.'
                % ('%.6f' % rho if rho is not None else 'unavailable'),
                RuntimeWarning, stacklevel=2)
        if not self.converged:
            ## ⚠ THIS USED TO BE SILENT.  `fsolve` builds the "No
            ## convergence" message and then discards it whenever
            ## `full_output=False`, which is how this call was written -- so
            ## a shooting solve that never converged returned a
            ## plausible-looking waveform with no diagnostic at all.  It was
            ## non-convergent on EVERY circuit, including a linear RLC whose
            ## answer was visibly close, which is why nobody noticed.
            warnings.warn(
                'PSS: the shooting solve did not converge in %d iterations '
                '(method=%r, %s). The returned waveform is the last iterate, '
                'not a periodic steady state -- raise maxiterations, or use '
                "method='euler', which solves a true Newton system while "
                'trapezoidal still uses successive substitution.'
                % (maxiterations, method,
                   'true Newton' if newton else 'successive substitution'),
                RuntimeWarning, stacklevel=2)
        
        X = [x0_ss]
        iq_last = None
        self._begin_period(x0_ss)
        for t in times:
            x = self.solve_timestep(X[-1], t, dt, iq_last=iq_last)
            iq_last = self._iq
            X.append(copy(x))

        X = toolkit.array(X[1:]).T

        # Insert reference node voltage
        X = toolkit.concatenate((X[:irefnode], 
                                 toolkit.zeros((1,len(times))), 
                                 X[irefnode:]))

        tpss = analysis.CircuitResult(self.cir, x=X, xdot=None,
                                      sweep_values=times, sweep_label='time', 
                                      sweep_unit='s')

        freqs, FX = freq_analysis(X[:,:-1], times[:-1])
        
        fpss = analysis.CircuitResult(self.cir, x=FX, xdot=None,
                                      sweep_values=freqs, sweep_label='freq', 
                                      sweep_unit='Hz')
        
        return InternalResultDict({'tpss': tpss, 'fpss': fpss})

class PAC(Analysis):
    """Small-signal analysis over a time varying operating point"""

    parameters  = [Parameter(name='analysis', desc='Analysis name', 
                             default='ac')]

    def __init__(self, cir, toolkit=None, **kvargs):
        self.parameters = super(PAC, self).parameters + self.parameters            
        super(PAC, self).__init__(cir, toolkit=toolkit, **kvargs)
    
    def solve(self, pss, freqs, refnode=gnd, period=1e-3, x0=None, timestep=1e-6,
              maxiterations=20):
        raise NotImplementedError(
            "PAC is withdrawn as unimplemented (stage 11). It forms the whole "
            "(N*M)x(N*M) operator densely: at N=137 unknowns and M=1000 time "
            "points that is 419.5 GiB (279.7 GiB for L in complex128 plus 139.8 "
            "GiB for B), which is not a tuning problem but a consequence of the "
            "formulation. It has also never been validated -- its only test was "
            "@unittest.skip('Skip failing test'). A thin advertised feature is "
            "worse than an absent one, so it says so instead of allocating.\n\n"
            "The body below is kept, unreachable, because it is the starting "
            "point for a matrix-free rewrite (Telichevesky, Kundert & White, DAC "
            "1995): PAC needs only products with the operator, never the operator "
            "itself. Use PSS for periodic steady state in the meantime.")

        tk = self.toolkit
        analysis_name = self.par.analysis
        print('solve PAC analysis_name = ' + analysis_name)
        ## Create U vector which is the RHS evaluated at every time instant
        T = pss.period
        times = pss.times[:-1]
        hs = tk.diff(pss.times)

        N = self.cir.n - 1 ## ref node removed
        M = len(times)

        irefnode = self.cir.get_node_index(refnode)
        (u0,) = remove_row_col((self.cir.u(0, analysis_name),), irefnode, self.toolkit)

        ## Create LHS matrix using backward Euler discretization
        L = tk.zeros((N*M, N*M),dtype=tk.cdouble)
        ## 0.1c: no dtype, so B was float64 while L is complex -- it worked by
        ## promotion, not by intent.  Corrected even though unreachable, so the
        ## starting point for a rewrite is not itself wrong.
        B = tk.zeros(L.shape, dtype=tk.cdouble)
        for i, (t, h, J, C) in enumerate(zip(times, hs, pss.Jtvec, pss.Cvec)):
            L[i*N:(i+1)*N, i*N:(i+1)*N] = J
            if i > 0:
                L[i*N:(i+1)*N, (i-1)*N:i*N] = -C / h
        B[0:N,(M-1)*N:M*N] = -pss.Cvec[-1] / hs[0]

        outfreq = []
        outV = []
        for fs in freqs:
            phase_shift = tk.zeros(N * M, dtype=tk.cdouble)
            u = tk.zeros(N * M, dtype=tk.cdouble)
            for i,t in enumerate(times):
                phase_shift[i*N:(i+1)*N] = tk.exp(2j*tk.pi*fs*t)
                u[i*N:(i+1)*N] = u0
            
            u *= phase_shift
            
            alpha = tk.exp(-2j*tk.pi*fs*T)

            ## Solve discrete-time AC-voltage vector
            v = tk.linearsolver(L + alpha*B, -u)
            
            ## multiply v matrix by exp(-j*2*pi*fs) so the spectrum
            ## is evaluated at 2*pi*(fs + 1/T) instead of 2*pi/T
            ## this will also make v T-periodic
            v_shifted = (v / phase_shift)

            freqs, V = freq_analysis(v_shifted.reshape(M,N),
                                     times, axis=0)

            outfreq.extend((abs(freqs + fs)).tolist())
            outV.extend(V.tolist())
            
        ## Sort on frequency
        freqs, X = zip(*sorted(zip(outfreq, outV)))

        X = np.array(X)
        freqs = np.array(freqs)

        # Insert reference node voltage
        irefnode = self.cir.get_node_index(refnode)
        X = tk.concatenate((X[:,:irefnode], 
                            tk.zeros((len(freqs),1)), 
                            X[:,irefnode:]), axis=1)


        res = analysis.CircuitResult(self.cir, x = X.T, 
                                        xdot=None,
                                        sweep_values=freqs, 
                                        sweep_label='freq', 
                                        sweep_unit='Hz')

        
        return res

        
            
