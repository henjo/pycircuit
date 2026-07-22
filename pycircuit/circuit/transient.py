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
        
        def limiter_func(xr, x0r):
            n = self.cir.n
            x = self.toolkit.zeros(n)
            x0_full = self.toolkit.zeros(n)
            
            x[:self.irefnode] = xr[:self.irefnode]
            x[self.irefnode+1:] = xr[self.irefnode:]
            
            x0_full[:self.irefnode] = x0r[:self.irefnode]
            x0_full[self.irefnode+1:] = x0r[self.irefnode:]
            
            self.cir.limit(x, x0_full, self.epar)
            return xr

        try:
            result = fsolve(refnode_removed(func, self.irefnode, self.toolkit), 
                            x0, 
                            full_output = True, 
                            reltol = self.par.reltol,
                            abstol = abstol, xtol=xtol,
                            maxiter = self.par.maxiter,
                            toolkit = self.toolkit,
                            limiter = limiter_func)
        except self.toolkit.linalg.LinAlgError as e:
            raise SingularMatrix(str(e))
        
        x, infodict, ier, mesg = result
        
        if ier != 1:
            raise NoConvergenceError(mesg)
        
        # Insert reference node voltage
        return self.toolkit.concatenate((x[:self.irefnode], self.toolkit.array([0.0]), x[self.irefnode:]))
    

    
    def get_diff(self, q, C, method=None):
        """Method used to calculate time derivative for charge storing elements (i_eq and g_eq)."""
        dt = self._dt
        method = method or self.par.method
        
        resultEuler = (q - self._qlast[0]) / dt
        
        # Bizzarri & Brambilla: Protect against sudden time step variations.
        # If the step shrinks by more than 10x, high-order history polynomials become invalid.
        # We forcefully drop order to Backward Euler to safely cross the discontinuity.
        if method == 'gear2' and not getattr(self, '_is_first_step', False):
            h1 = getattr(self, '_dt_last', dt)
            if dt / h1 < 0.1:
                method = 'euler'
        
        if getattr(self, '_is_first_step', False): #first step always requires backward euler
            geq = C / dt
            iq = resultEuler
            self._diff_error = self.toolkit.zeros(len(q))
        else:
            _, _, b_ = self._method.get('trap', (None, None, 0.5))
            resultTrap = 2*(q - self._qlast[0])/dt - self._iqlast[0]
            self._diff_error = resultTrap - resultEuler 
            if method == 'euler':
                iq = resultEuler
                geq = C / dt
            elif method in ('trapezoidal', 'trap'):
                iq = resultTrap
                geq = C / dt / b_
            elif method == 'gear2':
                h = dt
                h1 = getattr(self, '_dt_last', dt)
                alpha0 = (2*h + h1) / (h * (h + h1))
                alpha1 = -(h + h1) / (h * h1)
                alpha2 = h / (h1 * (h + h1))
                geq = C * alpha0
                iq = alpha0 * q + alpha1 * self._qlast[0] + alpha2 * self._qlast[1]
            else:
                a, b, b_ = self._method[method] 
                iq = (q - self.toolkit.dot(a, self._qlast[:len(a)])) / dt / b_ - self.toolkit.dot(b, self._iqlast[:len(b)]) / b_
                geq = C / dt / b_
                
        self._iq = iq 
        self._effective_method = method
        return iq, geq

    
    
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
    
    
    def solve(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6, provided_function=None, fixed_timestep=False, coupled_lte=False, analytical_eh=True):
        if coupled_lte:
            return self._solve_coupled(refnode, tend, x0, timestep, provided_function, analytical_eh)
        X = []
        self.irefnode=self.cir.get_node_index(refnode)
        n = self.cir.n
        if x0 is None:
            # Use DC operating point if x0 is not provided
            from pycircuit.circuit.dcanalysis import DC
            dc = DC(self.cir)
            try:
                dc_res = dc.solve()
                x0 = dc_res.x
            except Exception:
                x0 = self.toolkit.zeros(self.cir.n)
            x = x0
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
                    
                    if self._effective_method == 'euler':
                        Eg = -0.5 * (gn - gn_1)
                    elif self._effective_method in ('trap', 'trapezoidal'):
                        h = dt
                        h1 = getattr(self, '_dt_last', dt)
                        dd1 = (gn - gn_1) / h
                        dd2 = (gn_1 - gn_2) / h1
                        Eg = -(1.0/3.0) * h**2 * (dd1 - dd2) / (h + h1)
                    elif self._effective_method == 'gear2':
                        h = dt
                        h1 = getattr(self, '_dt_last', dt)
                        dd1 = (gn - gn_1) / h
                        dd2 = (gn_1 - gn_2) / h1
                        Eg = -(2.0/3.0) * h**2 * (dd1 - dd2) / (h + h1)
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
            self._dt_last = dt
            
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


    def _solve_coupled(self, refnode=gnd, tend=1e-3, x0=None, timestep=1e-6, provided_function=None, analytical_eh=True):
        import numpy as np
        from copy import copy
        from pycircuit.circuit.analysis import CircuitResult
        
        X = []
        self.irefnode = self.cir.get_node_index(refnode)
        n = self.cir.n
        if x0 is None:
            from pycircuit.circuit.dcanalysis import DC
            dc = DC(self.cir)
            try:
                dc_res = dc.solve()
                x0 = dc_res.x
            except Exception:
                x0 = self.toolkit.zeros(n)
        
        x = x0
        a,b,b_ = self._method[self.par.method]
        q0 = self.cir.q(x)
        self._qlast = self.toolkit.array([q0 for _ in range(max(2, len(a)))])
        self._iqlast = self.toolkit.zeros((max(2, len(b)), n))
        
        X.append(copy(x))
        timelist = []
        
        self._is_first_step = True
        t = 0.0
        h = timestep
        TRTOL = 7.0
        
        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))
        abstol = self.toolkit.concatenate((self.par.iabstol * ones_nodes,
                                         self.par.vabstol * ones_branches))
        reltol = self.par.reltol
        xtol = self.par.vabstol
        
        while t < tend:
            if t + h > tend:
                h = tend - t
                
            x_curr = copy(X[-1])
            h_curr = h
            converged = False
            
            for iter_idx in range(self.par.maxiter):
                # Compute F and J
                self._dt = h_curr
                C = self.cir.C(x_curr)
                q = self.cir.q(x_curr)
                iq, Geq = self.get_diff(q, C)
                F = self.cir.i(x_curr) + iq + self.cir.u(t + h_curr, analysis=self.par.analysis)
                J = self.cir.G(x_curr) + Geq
                
                F_r = self.toolkit.delete(F, self.irefnode)
                J_r = self.toolkit.delete(J, self.irefnode, axis=0)
                J_r = self.toolkit.delete(J_r, self.irefnode, axis=1)
                
                # Finite difference J_h
                eps = max(1e-8 * h_curr, 1e-15)
                self._dt = h_curr + eps
                iq_p, Geq_p = self.get_diff(q, C)
                F_p = self.cir.i(x_curr) + iq_p + self.cir.u(t + h_curr + eps, analysis=self.par.analysis)
                
                self._dt = h_curr - eps
                iq_m, Geq_m = self.get_diff(q, C)
                F_m = self.cir.i(x_curr) + iq_m + self.cir.u(t + h_curr - eps, analysis=self.par.analysis)
                
                J_h = (F_p - F_m) / (2 * eps)
                J_h_r = self.toolkit.delete(J_h, self.irefnode)
                self._dt = h_curr
                
                # Compute E, E_x, E_h
                def calc_E(x_val, h_val):
                    if self._is_first_step:
                        return 0.0
                    q_val = self.cir.q(x_val)
                    if self._effective_method in ("trapezoidal", "trap"):
                        dd2 = (q_val - self._qlast[0]) / h_val - self._iqlast[0]
                        dd2 = dd2 * 2.0 / h_val
                        lte = 1.0 / 12.0 * (h_val**3) * dd2
                    elif self._effective_method == "gear2":
                        dd1_n = (q_val - self._qlast[0]) / h_val
                        dd1_nm1 = (self._qlast[0] - self._qlast[1]) / self._dt_last
                        dd2_n = (dd1_n - dd1_nm1) / (h_val + self._dt_last)
                        lte = (h_val**2) * (h_val + self._dt_last) / 3.0 * dd2_n
                    else:
                        lte = self.toolkit.zeros(n)
                        
                    etol = reltol * self.toolkit.maximum(np.abs(q_val), np.abs(self._qlast[0])) + abstol
                    return np.max(np.abs(lte) / etol) - TRTOL
                
                E = calc_E(x_curr, h_curr)
                
                gamma_min = 0.7
                gamma_max = 3.0
                
                if analytical_eh and ((gamma_min - 1.0) * TRTOL <= E <= (gamma_max - 1.0) * TRTOL) and not self._is_first_step:
                    E_x_r = self.toolkit.zeros(len(F_r))
                    E_h = 1.0
                    E = 0.0  # Forces dh = 0
                elif self._is_first_step:
                    E_x_r = self.toolkit.zeros(len(F_r))
                    E_h = 1.0
                else:
                    if analytical_eh:
                        p = 3.0 if self._effective_method in ("trapezoidal", "trap") else 2.0
                        E_h = p * (E + TRTOL) / h_curr
                    else:
                        E_p = calc_E(x_curr, h_curr + eps)
                        E_m = calc_E(x_curr, h_curr - eps)
                        E_h = (E_p - E_m) / (2 * eps)
                    
                    E_x_r = self.toolkit.zeros(len(F_r))
                    if not analytical_eh:
                        for i in range(len(F_r)):
                            idx = i if i < self.irefnode else i + 1
                            x_p = copy(x_curr)
                            x_p[idx] += eps
                            E_xp = calc_E(x_p, h_curr)
                            
                            x_m = copy(x_curr)
                            x_m[idx] -= eps
                            E_xm = calc_E(x_m, h_curr)
                            
                            E_x_r[i] = (E_xp - E_xm) / (2 * eps)
                
                # Schur complement solve
                rhs = np.column_stack([-F_r, -J_h_r])
                try:
                    dx_res = self.toolkit.linearsolver(J_r, rhs)
                except Exception:
                    # Fallback if singular
                    dx_res = np.zeros_like(rhs)
                
                dx_0 = dx_res[:, 0]
                dx_h = dx_res[:, 1]
                
                if self._is_first_step:
                    dh = 0.0
                else:
                    denom = np.dot(E_x_r, dx_h) + E_h
                    if abs(denom) < 1e-20:
                        dh = 0.0
                    else:
                        dh = (-E - np.dot(E_x_r, dx_0)) / denom
                
                dh = max(-0.5 * h_curr, min(2.0 * h_curr, dh))
                
                dx_r = dx_0 + dx_h * dh
                
                x_next = copy(x_curr)
                x_next[:self.irefnode] += dx_r[:self.irefnode]
                x_next[self.irefnode+1:] += dx_r[self.irefnode:]
                
                self.cir.limit(x_next, x_curr, self.epar)
                
                x_curr = x_next
                h_curr += dh
                
                x_curr_r = self.toolkit.delete(x_curr, self.irefnode)
                abstol_r = self.toolkit.delete(abstol, self.irefnode)
                
                converged_x = np.all(np.abs(dx_r) < reltol * np.maximum(np.abs(x_curr_r), 1e-12) + abstol_r)
                converged_h = np.abs(dh) < 0.15 * h_curr
                
                # NEW KCL Check for Option A
                I_scale_r = np.dot(np.abs(J_r), np.abs(x_curr_r)) + np.abs(F_r)
                converged_F = np.all(np.abs(F_r) < reltol * I_scale_r + abstol_r)
                
                if converged_x and converged_h and converged_F:
                    converged = True
                    break
                    
            if not converged:
                # If NR failed, reduce step drastically
                h *= 0.25
                if h < 1e-15:
                    raise RuntimeError("Timestep too small")
                continue
                
            t += h_curr
            h = h_curr
            timelist.append(t)
            X.append(copy(x_curr))
            
            self._dt = h_curr
            self._dt_last = h_curr
            self._is_first_step = False
            self._iqlast = self.toolkit.concatenate((self.toolkit.array([iq]), self._iqlast))[:-1]
            self._qlast = self.toolkit.concatenate((self.toolkit.array([q]), self._qlast))[:-1]
            
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
