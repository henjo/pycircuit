# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.circuit.analysis import *
from pycircuit.circuit.dcanalysis import DC
import theanotk

class Transient(Analysis):
    """Simple transient analysis class.

    Time step is fixed.

    i(t) = c*dv/dt
    v(t) = L*di/dt

    The usual companion models are used.
    backward euler:
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
    >>> expected = 6.3
    >>> abs(res.v(n2, gnd)[-1] - expected) < 1e-2*expected #node 2 of last x
    True

    Linear circuit example:
    >>> from pycircuit.circuit.elements import ISin
    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> c['Isin'] = ISin(gnd, n1, ia=1e-3, freq=16e3)
    >>> c['R'] = R(n1, gnd, r=200)
    >>> c['C'] = C(n1, gnd, c=1e-6)
    >>> c['L'] = L(n1, gnd, L=1e-4)
    >>> tran = Transient(c)
    >>> res = tran.solve(tend=260e-6,timestep=1e-6)
    >>> expected = 0.063
    >>> abs(res.v(n1,gnd)[-1]) < 1e-1*expected #node 2 of last x
    True
    
    """
    ## TODO:
    ## * Implement automatic timestep adjustment, using difference between
    ##   BE and trapezoidal as a measure of the error.
    ##   Reference: "Time Step Control in Transient Analysis", by SHUBHA VIJAYCHAND
    

    parameters = Analysis.parameters + \
        [Parameter(name='analysis', desc='Analysis name', 
                   #default='transient'),
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

    def __init__(self, cir, toolkit=theanotk, **kvargs):
        self.parameters = super(Transient, self).parameters + self.parameters            
        super(Transient, self).__init__(cir, toolkit=toolkit, **kvargs)
    
        self._method={
            "euler":(self.toolkit.array([1.]),self.toolkit.array([0.]),1.),
            "trap":(self.toolkit.array([1.]),self.toolkit.array([0.5]),0.5),
            "trapezoidal":(self.toolkit.array([1.]),self.toolkit.array([0.5]),0.5),
            "gear2":(self.toolkit.array([4./3,-1./3]),self.toolkit.array([0]),2./3)
            }
        self._qlast  = None #q history
        self._iqlast = None #dq/dt history
        
        self._dt = None
        self._diff_error = None #used for saving difference between euler and trapezoidal

        self.epar.analysis = 'tran'
    
    ## This is borrowed from dcanalysis.py, would like to 
    ## import it from there instead.
    ## But it's an object method requiring a DC as self
    ## so using DC._newton doesn't work
    def _newton(self, func, x0): 
        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))
        
        abstol = self.toolkit.concatenate((self.par.iabstol * ones_nodes,
                                 self.par.vabstol * ones_branches))
        xtol = self.toolkit.concatenate((self.par.vabstol * ones_nodes,
                                 self.par.iabstol * ones_branches))
        
        try:
            result = fsolve(func, 
                            x0, 
                            full_output = True, 
                            reltol = self.par.reltol,
                            abstol = abstol, xtol=xtol,
                            maxiter = self.par.maxiter,
                            toolkit = self.toolkit)
        except self.toolkit.linalg.LinAlgError, e:
            raise SingularMatrix(e.message)
        
        x, infodict, ier, mesg = result
        
        if ier != 1:
            raise NoConvergenceError(mesg)
        
        # Insert reference node voltage
        return x
    
    def get_timestep(self,endtime,dtmin=1e-12):
        """Method to provide the next timestep for transient simulation.
        
        """
        ## Use the difference between euler and trapezoidal
        ## Error as abs(self._diff_error)/abs(self._iq) - iq_tolerance
        ## If error > 0: decrease timestep
        ## If error < 0: increase timestep
        
        ## Start with simple variant:
        ## increase as dt *=2
        ## decrease as dt /=2
        ## with convergence error, both decrease and reject timestep
        ## support for rejecting timesteps need to be implemented
        
        iq_tolerance = 1e-2
        iq_p = 1e4 #proportionality constant
        iq_i = 1 #integrator constant
        dt=self._dt
        t=0
        while t<endtime:
            yield t,dt
            de=self._diff_error
            iq=self._iq
            if (de != None) and (iq != None):
                #iq_error=self.toolkit.dot(de,de)/self.toolkit.dot(iq,iq)-iq_tolerance
                #print iq_error
                dt = max(dt, dtmin)
            t+=dt
    
    def get_diff(self,q,C):#shouldn't I provide an x0 here?
        """Method used to calculate time derivative for charge storing elements (i_eq and g_eq).
        
        Calculates approximate derivatives, both for backward euler and trapezoidal. 
        The difference between these can be used to determine the next timestep (or 
        reject the last). The difference is stored in a class variable/attribute and
        return value is one of the calculated derivatives, dependent on selected
        integration method.
        """
        #calculate in a more general way with coefficients dependent on method
        #the amount of history values is determined by the length of the coefficient-vector
        
        dt=self._dt
        a,b,b_=self._method[self.par.method] 
        resultEuler = (q-self._qlast[0])/dt
        if self._iqlast == None: #first step always requires backward euler
            geq=C/dt
            n=self.cir.n
            self._iqlast=self.toolkit.zeros((len(b),n)) #initialize history vectors at first step
            iq = resultEuler
        else:
            geq=C/dt/b_
            resultTrap = 2*(q-self._qlast[0])/dt-self._iqlast[0]
            self._diff_error = resultTrap-resultEuler # Difference between euler and trap.
            if self.par.method == 'euler':
                iq = resultEuler
            elif self.par.method == 'trapezoidal':
                iq = resultTrap
            else:
                iq=(q-self.toolkit.dot(a,self._qlast))/dt/b_ - self.toolkit.dot(b,self._iqlast)/b_
        self._iq=iq #make accessible by get_timestep
        return iq,geq
    
    
    def solve_timestep(self, x0, t, refnode=gnd, provided_function=None):
        #if provided_function is not None, it is called as a function with 
        #most of what is calculated during a time_step, f,J,ueq,Geq,xlast,x
        
        self.epar.t = t
        
        n=self.cir.n
        x0 = x0
        dt = self._dt
        
        def func(x):
            self.na.update(x, self.epar)
            C = self.na.C
            q = self.na.q
            iq,Geq = self.get_diff(q,C)
            f = self.na.i + iq + self.na.u
            J = self.na.G + Geq #return C somehow?
            return self.toolkit.array(f, dtype=float), self.toolkit.array(J, dtype=float)
        
        x=self._newton(func,x0)
        #history update
        self._iqlast = self.toolkit.concatenate((self.toolkit.array([self._iq]),self._iqlast))[:-1]
        self._qlast = self.toolkit.concatenate((self.toolkit.array([self.na.q]),self._qlast))[:-1]
        
        # Insert reference node voltage
        #x = self.toolkit.concatenate((x[:irefnode], self.toolkit.array([0.0]), x[irefnode:]))
        if provided_function != None:
            result=x,provided_function(f,J,C)
        else:
            result=x,None

        return result
    
    
    def solve(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6, provided_function=None):
        #provided_function is a function that is sent to solve_timestep for evaluation
        self.na.setup()
        
        X = [] # will contain a list of all x-vectors
        n = self.na.n
        self._dt = timestep
        if x0 is None:
            x = self.toolkit.zeros(n)
        else:
            x = x0 
        
        a,b,b_=self._method[self.par.method] 
        self._qlast=self.toolkit.zeros((len(a),n))#initialize q-history vector
        #shift in q(x0) to q-history
        self._qlast = self.toolkit.concatenate((self.toolkit.array([self.na.q]),self._qlast))[:-1]
        #is this still needed
        order=1 #number of past x-values needed
        for i in xrange(order):
            X.append(copy(x))
        
        times = self.get_timestep(tend)
        timelist=[] #for plotting purposes
        self._iqlast=None #forces first step to be Backward Euler
        for t,dt in times:
            timelist.append(t)
            self._dt=dt
            x,feval=self.solve_timestep(X[-1], t, provided_function=provided_function)
            X.append(copy(x))
        X = self.toolkit.array(X[1:]).T
        timelist = self.toolkit.array(timelist)
        
        #print("steps: "+str( len(timelist)))
        
        self.result = CircuitResult(self.na, x=X, xdot=None,
                                    sweep_values=timelist, 
                                    sweep_label='time', 
                                    sweep_unit='s')
        
        return self.result


if __name__ == "__main__":
    import doctest
    doctest.testmod()
