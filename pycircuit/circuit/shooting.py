from pycircuit.post import InternalResultDict
from circuit import gnd
from analysis import Analysis, remove_row_col
from copy import copy
import analysis
import numpy as np

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
            ueq = -toolkit.dot(Geq,xlast)
            f =  self.cir.i(x) + self.cir.q(x)/dt + self.cir.u(t) + ueq
            J = self.cir.G(x) + Geq
            (f,J,C) = remove_row_col((f,J,C), irefnode)
            self._Jf, self._C = J, C
            return f, J

        rtol = 1e-4

        x = analysis.fsolve(func, x0, reltol=rtol)
        # Insert reference node voltage
        #x = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
        return x


    def solve(self, refnode=gnd, period=1e-3, x0=None, timestep=1e-6, 
              maxiterations=20):
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
            for t in times[1:]:
                x = copy(self.solve_timestep(x, t, dt))
                Jshoot = np.mat(self._Jf).I * C * Jshoot
                C = copy(np.mat(self._C))

            residual = x0 - x

            D = np.mat(toolkit.eye(n-1))
            return residual, D - alpha * Jshoot
        
        ## Find periodic steady state x-vector
        x0_ss = analysis.fsolve(func, x, maxiter=maxiterations)
        
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

        npoints = len(times) - 1
        if X.dtype is np.complex:
            FX = np.fft.fftshift(np.fft.fft(X[:,:-1], axis=-1)) / npoints
            freqs = np.fft.fftshift(np.fft.fftfreq(npoints, d=dt))
        else:
            freqs = np.fft.fftfreq(npoints, d=dt)[:np.ceil(npoints / 2.)]
            FX = np.fft.fft(X[:,:-1], axis=-1)[:,:len(freqs)] / npoints
            ## Fold energy from negative frequencies
            FX[:,1:] *= np.sqrt(2)

        fpss = analysis.CircuitResult(self.cir, x=FX, xdot=None,
                                      sweep_values=freqs, sweep_label='freq', 
                                      sweep_unit='Hz')
        
        return InternalResultDict({'tpss': tpss, 'fpss': fpss})
