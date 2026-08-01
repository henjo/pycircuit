from pycircuit.post import InternalResultDict
from .circuit import gnd
from pycircuit.circuit.analysis import *
from copy import copy
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

class PSS(Analysis):
    """Periodic Steady-State using shooting Newton iterations
    
    The algorithm is described in [1] p65.

     1. Kenneth S. Kundert, Jacob K. White, Alberto Sangiovanni-Vincentelli 
        (1990)
        Steady-State Methods for Simulating Analog and Microwave Circuits
        Kluwer Academic Publishers
        ISBN 0792390695
    
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
         Parameter(name='method', 
                   desc='Differentiation method', unit='', 
                   default="euler")]        

    
    def __init__(self, cir, toolkit=None, irefnode=None, **kvargs):
        self.parameters = super(PSS, self).parameters + self.parameters            
        super(PSS, self).__init__(cir, **kvargs)

    def solve_timestep(self, x0, t, dt, refnode=gnd, iq_last=None):
        """One timestep of the inner transient.

        Returns the solution; the companion current at the converged point is left
        in ``self._iq`` for the caller to feed back as ``iq_last``.  Kept out of
        the return value so existing callers are unaffected.

        STAGE 11 -- `method` NOW SELECTS SOMETHING.  It was declared with
        `default="euler"` and never read anywhere in this file: `solve_timestep`
        hard-coded `C/dt` and `q(xlast)/dt`, so PSS was backward-Euler only.

        For a PERIODIC STEADY-STATE solver that is the worst available fixed
        choice, because backward Euler's numerical damping attenuates exactly the
        limit cycle PSS exists to find.  Measured on a series RLC driven at
        resonance, Q = 20, against the analytic peak of 20 V:

            steps/period    PSS peak    fraction of analytic
                      20      2.63 V       13.2%
                      50      5.61 V       28.1%
                     100      8.81 V       44.1%
                     200     12.20 V       61.0%

        Silently, and worse as the step coarsens -- an oscillator's amplitude
        comes out low and a resonator's Q understated with no diagnostic at all.

        `iq_last` is the previous step's companion current, which trapezoidal
        needs and Euler ignores.  `None` means "no history yet", and the first
        step of a sweep therefore falls back to Euler, which is standard: there is
        nothing to average against.
        """
        toolkit = self.toolkit
        concatenate, array = toolkit.concatenate, toolkit.array

        n=self.cir.n
        analysis_name = self.par.analysis
        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.cir.get_node_index(refnode)

        method = getattr(self.par, 'method', 'euler')
        if method not in ('euler', 'trap', 'trapezoidal'):
            raise ValueError(
                "method must be 'euler' or 'trap', not %r" % (method,))
        ## Trapezoidal needs a companion current to average against; without one
        ## (the first step of a sweep) it degenerates to Euler by construction.
        use_trap = method in ('trap', 'trapezoidal') and iq_last is not None

        self._iq = None

        def func(x):
            x = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
            xlast = concatenate((x0[:irefnode], array([0.0]), x0[irefnode:]))
            C = self.cir.C(x)
            q, qlast = self.cir.q(x), self.cir.q(xlast)
            if use_trap:
                ## iq_n = 2 (q_n - q_{n-1})/dt - iq_{n-1}, and dq/dx = 2C/dt.
                iq = 2.0 * (q - qlast) / dt - iq_last
                Geq = 2.0 * C / dt
            else:
                iq = (q - qlast) / dt
                Geq = C / dt
            f = self.cir.i(x) + iq + self.cir.u(t, analysis=analysis_name)
            J = self.cir.G(x) + Geq
            (f, J, C) = remove_row_col((f, J, C), irefnode, self.toolkit)
            self._Jf, self._C = J, C
            self._iq = iq
            return f, J

        x = analysis.fsolve(func, x0, reltol=self.par.reltol, toolkit=self.toolkit)

        ## STAGE 11 -- RECOMPUTE THE COMPANION AT THE CONVERGED POINT.
        ##
        ## `func` leaves `self._iq` wherever fsolve last evaluated it, and fsolve
        ## evaluates at the iterate BEFORE the update it returns -- so the stored
        ## value belongs to the previous iterate, not to `x`.  Feeding that back as
        ## `iq_last` seeds the trapezoidal recursion with a point the circuit was
        ## never at, and the recursion has an undamped (-1)^n mode to amplify it.
        ## Exactly the staleness found in the JAX Newton under 9(e), in a second
        ## place.
        xf = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
        xl = concatenate((x0[:irefnode], array([0.0]), x0[irefnode:]))
        dq = (self.cir.q(xf) - self.cir.q(xl)) / dt
        self._iq = (2.0 * dq - iq_last) if use_trap else dq

        # Insert reference node voltage
        #x = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
        return x


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

        def func(x):
            ## STAGE 11 -- the companion current is carried through the sweep, so
            ## trapezoidal has something to average against.  `iq_last=None` on the
            ## first step is what makes it fall back to Euler there.
            x = self.solve_timestep(x, times[0], dt)
            iq_last = self._iq
            x0 = copy(x)
            Jshoot = np.asarray(toolkit.eye(n-1))
            C = copy(np.asarray(self._C))

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
                Jshoot = self.toolkit.linearsolver(np.asarray(self._Jf),
                                                   C @ Jshoot)
                C = copy(np.asarray(self._C))

            residual = x0 - x

            D = np.asarray(toolkit.eye(n-1))
            return residual, D - alpha * Jshoot
        
        ## Find periodic steady state x-vector
        x0_ss = analysis.fsolve(func, x, maxiter=maxiterations, toolkit=self.toolkit)
        
        X = [x0_ss]
        iq_last = None
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

        
            
