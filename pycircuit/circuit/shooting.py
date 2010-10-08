from pycircuit.post import InternalResultDict
from circuit import gnd
from analysis import Analysis, remove_row_col
from copy import copy
import analysis
import numpy as np

def freq_analysis(x, t, rms = True, axis=-1, freqoffset = 0):
    """Return dft of equidistant sampled signal x"""
    
    npoints = np.size(x, axis)

    dt = t[1] - t[0]

    if x.dtype in (np.complex128, np.complex):
        X = np.fft.fftshift(np.fft.fft(x, axis=axis),axes=(axis,)) / npoints
        freqs = np.fft.fftshift(np.fft.fftfreq(npoints, d=dt))
    else:
        freqs = np.fft.fftfreq(npoints, d=dt)[:np.ceil(npoints / 2.)]
        slices = [slice(None)] * x.ndim
        slices[axis] = slice(0, len(freqs))
        X = np.fft.fft(x, axis=axis)[slices] / npoints
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
    def solve_timestep(self, x0, t, dt, refnode=gnd):
        toolkit = self.toolkit
        concatenate, array = toolkit.concatenate, toolkit.array

        n=self.cir.n

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.cir.get_node_index(refnode)

        def func(x):
            x = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
            xlast = concatenate((x0[:irefnode], array([0.0]), x0[irefnode:]))
            C = self.cir.C(x)
            Geq = C/dt
            ueq = -self.cir.q(xlast)/dt
            f =  self.cir.i(x) + self.cir.q(x)/dt + self.cir.u(t) + ueq
            J = self.cir.G(x) + Geq
            (f,J,C) = remove_row_col((f,J,C), irefnode, self.toolkit)
            self._Jf, self._C = J, C
            return f, J

        rtol = 1e-4

        x = analysis.fsolve(func, x0, reltol=rtol, toolkit=self.toolkit)
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
            x = self.solve_timestep(x, times[0], dt)
            x0 = copy(x)
            Jshoot = np.mat(toolkit.eye(n-1))
            C = copy(np.mat(self._C))

            ## Save C and transient jacobian for PAC analysis
            self.Cvec = [copy(self._C)]
            self.Jtvec = [copy(self._Jf)]
            self.times = times
            for t in times[1:]:
                x = copy(self.solve_timestep(x, t, dt))
                self.Cvec.append(copy(self._C))
                self.Jtvec.append(copy(self._Jf))
                Jshoot = np.mat(self._Jf).I * C * Jshoot
                C = copy(np.mat(self._C))

            residual = x0 - x

            D = np.mat(toolkit.eye(n-1))
            return residual, D - alpha * Jshoot
        
        ## Find periodic steady state x-vector
        x0_ss = analysis.fsolve(func, x, maxiter=maxiterations, toolkit=self.toolkit)
        
        X = [x0_ss]
        for t in times:
            x=self.solve_timestep(X[-1],t,dt)
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

    def solve(self, pss, freqs, refnode=gnd, period=1e-3, x0=None, timestep=1e-6, 
              maxiterations=20):
        tk = self.toolkit

        ## Create U vector which is the RHS evaluated at every time instant
        T = pss.period
        times = pss.times[:-1]
        hs = tk.diff(pss.times)

        N = self.cir.n - 1 ## ref node removed
        M = len(times)

        irefnode = self.cir.get_node_index(refnode)
        (u0,) = remove_row_col((self.cir.u(0, analysis='ac'),), irefnode, self.toolkit)

        ## Create LHS matrix using backward Euler discretization
        L = tk.zeros((N*M, N*M),dtype=tk.complex)
        B = tk.zeros(L.shape)
        for i, (t, h, J, C) in enumerate(zip(times, hs, pss.Jtvec, pss.Cvec)):
            L[i*N:(i+1)*N, i*N:(i+1)*N] = J
            if i > 0:
                L[i*N:(i+1)*N, (i-1)*N:i*N] = -C / h
        B[0:N,(M-1)*N:M*N] = -pss.Cvec[-1] / hs[0]

        outfreq = []
        outV = []
        for fs in freqs:
            phase_shift = tk.zeros(N * M, dtype=complex)
            u = tk.zeros(N * M, dtype=complex)
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

        
            
