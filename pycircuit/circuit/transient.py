# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.circuit.analysis import *

class Transient(Analysis):
    """Simple transient analysis class.

    Time step is fixed.
    Only method is currently backward euler.

    i(t) = c*dv/dt
    v(t) = L*di/dt

    backward euler:
    The usual companion models are used.
    i(n+1) = c/dt*(v(n+1) - v(n)) = geq*v(n+1) + Ieq
    v(n+1) = L/dt*(i(n+1) - i(n)) = req*i(n+1) + Veq

    def F(x): return i(x)+Geq(x)*x+u(x)+ueq(x0), G(x)+Geq(x)
    x0=x(n)
    x(n+1) = fsolve(F, x0, fprime=J)

    Linear circuit example:
    >>> circuit.default_toolkit = numeric
    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> n2 = c.add_node('net2')
    >>> c['Is'] = IS(gnd, n1, i=10)    
    >>> c['R1'] = R(n1, gnd, r=1)
    >>> c['R2'] = R(n1, n2, r=1e3)
    >>> c['R3'] = R(n2, gnd, r=100e3)
    >>> c['C'] = C(n2, gnd, c=1e-5)
    >>> tran = Transient(c)
    >>> res = tran.solve(tend=10e-3,timestep=1e-4)
    >>> res.v(n2, gnd)[-1] #node 2 of last x
    6.3

    Linear circuit example:
    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> c['Is'] = IS(gnd, n1, i=10e-3)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c['C'] = C(n1, gnd, c=1e-5)
    >>> c['L'] = L(n1, gnd, L=1e-3)
    >>> tran = Transient(c)
    >>> res = tran.solve(tend=150e-6,timestep=1e-6)
    >>> res.v(n1,gnd)[-1]
    0.99

    """
    ## TODO:
    ## * Make a timestep method that provides the next timestep
    ## * Implement trapezoidal method
    ## * Implement automatic timestep adjustment, using difference between
    ## BE and trapezoidal as a measure of the error.
    ## Reference: "Time Step Control in Transient Analysis", by SHUBHA VIJAYCHAND
    ## Perhaps use a PID-regulator for timestep adustment

    def get_timestep():
        """Method to provide the next timestep for transient simulation.
        
        """
        pass

    def get_diff(self,x,xlast,dt,iqlast=None):
        """Method used to calculate time derivative for charge storing elements.

        Calculates and returns approximate derivatives, both for backward euler 
        and trapezoidal. The difference between these can be used to determine
        the next timestep (or reject the last).
        """
        #BE: i(x)=(q(x)-q(xlast))/dt
        #Trap: i(x)=(q(x)-q(xlast))*2/dt-i_qlast
        q=self.cir.q(x)
        qlast=self.cir.q(xlast)
        #return 2*(qx-qxlast)/dt-self.cir.i_qlast
        return (q-qlast)/dt

    
    def solve_timestep(self, x0, t, dt, iqlast=None, refnode=gnd, rtol=1e-4, provided_function=None):
        #if provided_function is not None, it is called as a function with 
        #most of what is calculated during a time_step, f,J,ueq,Geq,xlast,x

        n=self.cir.n
        x0 = x0
        dt = dt

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.cir.get_node_index(refnode)

        def func(x):
            x = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
            xlast = concatenate((x0[:irefnode], array([0.0]), x0[irefnode:]))
            C = self.cir.C(x)
            Geq = C/dt
            ##Store dynamic current iq, so it can be reached in solve
            self._iq = self.get_diff(x,xlast,dt,iqlast=iqlast)
            f =self.cir.i(x) + self._iq + self.cir.u(t)
            J = self.cir.G(x) + Geq
            f, J, C = remove_row_col((f,J,C), irefnode)
            return array(f, dtype=float), array(J, dtype=float)

        x = fsolve(func, x0)
        # Insert reference node voltage
        #x = concatenate((x[:irefnode], array([0.0]), x[irefnode:]))
        if provided_function != None:
            result=x,provided_function(f,J,C)
        else:
            result=x,None
        return result


    def solve(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6, rtol=1e-4,provided_function=None):
        #should perhaps call with an option dictionary for things like rtol, method etc.
        #provided_function is a function that is sent to solve_timestep for evaluation

        X = [] # will contain a list of all x-vectors
        irefnode=self.cir.get_node_index(refnode)
        n = self.cir.n
        dt = timestep
        if x0 is None:
            x = zeros(n-1) #currently without reference node !
        else:
            x = x0 # reference node not included !
        order=1 #number of past x-values needed
        for i in xrange(order):
            X.append(copy(x))

        #create vector with timepoints and a more fitting dt
        times,dt=np.linspace(0,tend,num=int(tend/dt),endpoint=True,retstep=True)
        
        #initial step (t=0), always use backward euler
        x,feval = self.solve_timestep(X[-1],times[0],dt,rtol=rtol,provided_function=provided_function)
        X.append(copy(x))
        # Store last dynamic current (iq) for use in trapezoidal method
        # calculated in solve_timestep by method get_diff
        iqlast = self._iq

        for t in times[1:]: #all times but first (0)
            x,feval=self.solve_timestep(X[-1],t,dt,rtol=rtol,provided_function=provided_function,iqlast=iqlast)
            X.append(copy(x))
            #save iqlast for Trapezoidal method
            iqlast = self._iq
        X = self.toolkit.array(X[1:]).T

        # Insert reference node voltage
        X = self.toolkit.concatenate((X[:irefnode], 
                                      self.toolkit.zeros((1,len(times))), 
                                      X[irefnode:]))

        self.result = CircuitResult(self.cir, x=X, xdot=None,
                                    sweep_values=times, 
                                    sweep_label='time', 
                                    sweep_unit='s')
        
        return self.result


if __name__ == "__main__":
    import doctest
    doctest.testmod()