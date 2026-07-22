# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from pycircuit.circuit.analysis import *
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.dcanalysis import refnode_removed

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

    def __init__(self, cir, toolkit=None, irefnode=None, **kvargs):
        self.parameters = super(Transient, self).parameters + self.parameters            
        super(Transient, self).__init__(cir, **kvargs)
    
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
        
        (x0, abstol, xtol) = remove_row_col((x0, abstol, xtol), self.irefnode, self.toolkit)
        
        try:
            result = fsolve(refnode_removed(func, self.irefnode, self.toolkit), 
                            x0, 
                            full_output = True, 
                            reltol = self.par.reltol,
                            abstol = abstol, xtol=xtol,
                            maxiter = self.par.maxiter,
                            toolkit = self.toolkit)
        except self.toolkit.linalg.LinAlgError as e:
            raise SingularMatrix(e.message)
        
        x, infodict, ier, mesg = result
        
        if ier != 1:
            raise NoConvergenceError(mesg)
        
        # Insert reference node voltage
        return self.toolkit.concatenate((x[:self.irefnode], self.toolkit.array([0.0]), x[self.irefnode:]))
    

    
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
        
        if getattr(self, '_is_first_step', False): #first step always requires backward euler
            geq=C/dt
            iq = resultEuler
            self._diff_error = self.toolkit.zeros(len(q))
        else:
            geq=C/dt/b_
            resultTrap = 2*(q-self._qlast[0])/dt-self._iqlast[0]
            self._diff_error = resultTrap-resultEuler # Difference between euler and trap.
            if self.par.method == 'euler':
                iq = resultEuler
            elif self.par.method == 'trapezoidal':
                iq = resultTrap
            else:
                iq=(q-self.toolkit.dot(a,self._qlast[:len(a)]))/dt/b_ - self.toolkit.dot(b,self._iqlast[:len(b)])/b_
        self._iq=iq #make accessible by get_timestep
        return iq,geq
    
    
    def solve_timestep(self, x0, t, refnode=gnd, provided_function=None):
        n=self.cir.n
        dt = self._dt
        
        def func(x):
            C = self.cir.C(x)
            q=self.cir.q(x)
            iq,Geq = self.get_diff(q,C)
            f =self.cir.i(x) + iq + self.cir.u(t, analysis=self.par.analysis)
            J = self.cir.G(x) + Geq
            return self.toolkit.array(f, dtype=float), self.toolkit.array(J, dtype=float)
        
        x=self._newton(func,x0)
        f, J = func(x)
        
        if provided_function is not None:
            result=x,provided_function(f,J,self.cir.C(x)), J, f
        else:
            result=x,None, J, f
        return result
    
    
    def solve(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6, provided_function=None, fixed_timestep=False):
        X = []
        self.irefnode=self.cir.get_node_index(refnode)
        n = self.cir.n
        if x0 is None:
            # Calculate DC operating point as initial state
            dc = DC(self.cir, toolkit=self.toolkit, refnode=refnode)
            dc.par.reltol = self.par.reltol
            dc.par.vabstol = self.par.vabstol
            dc.par.iabstol = self.par.iabstol
            dc.par.maxiter = self.par.maxiter
            dc_res = dc.solve()
            x = dc_res.x
        else:
            x = x0 
        
        a,b,b_=self._method[self.par.method] 
        q0 = self.cir.q(x)
        self._qlast = self.toolkit.array([q0 for _ in range(max(2, len(a)))])
        self._iqlast = self.toolkit.zeros((max(2, len(b)), n))
        
        X.append(copy(x))
        
        timelist = []
        self._is_first_step = True
        t = 0.0
        max_step = timestep
        dt = timestep
        TRTOL = 7.0
        
        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))
        abstol = self.toolkit.concatenate((self.par.iabstol * ones_nodes,
                                         self.par.vabstol * ones_branches))

        while t < tend:
            if t + dt > tend:
                dt = tend - t
            
            self._dt = dt
            next_t = t + dt
            
            try:
                x, feval, J, f = self.solve_timestep(X[-1], next_t, provided_function=provided_function)
            except NoConvergenceError:
                dt = dt * 0.25
                continue
                
            if not fixed_timestep:
                # LTE Computation (Yao et al. ICECS 2014)
                if self._is_first_step:
                    err = 0.5 
                else:
                    gn = self._iq
                    gn_1 = self._iqlast[0]
                    gn_2 = self._iqlast[1]
                    
                    if self.par.method == 'euler':
                        Eg = -0.5 * (gn - gn_1)
                    elif self.par.method in ('trap', 'trapezoidal'):
                        Eg = -(1./6.) * (gn - 2*gn_1 + gn_2)
                    elif self.par.method == 'gear2':
                        Eg = -(1./3.) * (gn - 2*gn_1 + gn_2)
                    else:
                        Eg = self.toolkit.zeros(len(gn))
                    
                    J_reduced, Eg_reduced = remove_row_col((J, Eg), self.irefnode, self.toolkit)
                    
                    try:
                        lte_reduced = self.toolkit.linalg.solve(J_reduced, Eg_reduced)
                    except Exception:
                        lte_reduced = Eg_reduced
                    
                    lte = self.toolkit.concatenate((lte_reduced[:self.irefnode], self.toolkit.array([0.0]), lte_reduced[self.irefnode:]))
                    
                    import numpy as np
                    etol = self.par.reltol * np.maximum(np.abs(x), np.abs(X[-1])) + abstol
                    err_array = np.abs(lte) / etol
                    err = np.max(err_array)
                
                    if err > 1.0:
                        dt = dt * max(0.2, (1.0 / err)**0.5)
                        continue
            
            t = next_t
            timelist.append(t)
            X.append(copy(x))
            
            self._iqlast = self.toolkit.concatenate((self.toolkit.array([self._iq]), self._iqlast))[:-1]
            self._qlast = self.toolkit.concatenate((self.toolkit.array([self.cir.q(x)]), self._qlast))[:-1]
            
            self._is_first_step = False
            
            # Predict next step size
            if not fixed_timestep:
                dt = dt * min(2.0, (TRTOL / max(err, 1e-12))**0.5)
                dt = min(dt, max_step)
            
        X = self.toolkit.array(X[1:]).T
        timelist = self.toolkit.array(timelist)
        
        self.result = CircuitResult(self.cir, x=X, xdot=None,
                                    sweep_values=timelist, 
                                    sweep_label='time', 
                                    sweep_unit='s')
        
        return self.result


if __name__ == "__main__":
    import doctest
    doctest.testmod()
