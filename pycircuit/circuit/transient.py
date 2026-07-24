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
    def _get_integrator(self):
        from pycircuit.circuit.integrator import EulerIntegrator, TrapezoidalIntegrator, Gear2Integrator
        if self.par.method == 'euler':
            return EulerIntegrator()
        elif self.par.method in ('trap', 'trapezoidal'):
            return TrapezoidalIntegrator()
        elif self.par.method == 'gear2':
            return Gear2Integrator()
        else:
            raise ValueError(f"Unknown integration method: {self.par.method}")
    

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
                   default="euler"),
         Parameter(name='uic',
                   desc='Use initial conditions (skip DC OP computation)', unit='',
                   default=False),
         Parameter(name='minbreak',
                   desc='Minimum time difference for breakpoint events', unit='s',
                   default=1e-14),
         Parameter(name='bypass',
                   desc='Enable device model bypassing', unit='',
                   default=False),
         Parameter(name='minstep',
                   desc='Minimum timestep to prevent infinite loops', unit='s',
                   default=1e-18),
         Parameter(name='bypasstol',
                   desc='Bypass tolerance for device models', unit='V',
                   default=None)]

    def __init__(self, cir, toolkit=None, irefnode=None, **kvargs):
        self.parameters = super(Transient, self).parameters + self.parameters            
        super(Transient, self).__init__(cir, **kvargs)
    
        self._qlast  = None #q history
        self._iqlast = None #dq/dt history
        
        self._dt = None
        self._dt_last = None
        self._is_first_step = True
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
            x = self.toolkit.insert(xr, self.irefnode, 0.0)
            x0_full = self.toolkit.insert(x0r, self.irefnode, 0.0)
            
            x = self.cir.limit(x, x0_full, self.epar)
            return self.toolkit.concatenate((x[:self.irefnode], x[self.irefnode+1:]))

        from pycircuit.circuit.nrsolver import NoConvergenceError
        solver = self._get_nrsolver()
        scaler = self._get_scaler()
        try:
            x_res, _ = solver.solve_system(
                x0,
                refnode_removed(func, self.irefnode, self.toolkit),
                self.toolkit,
                self.par.reltol,
                abstol,
                xtol,
                self.par.maxiter,
                limiter=limiter_func,
                scaler=scaler
            )
        except Exception as e:
            if "Singular" in str(e) or "linearsolver" in str(e).lower() or "LinAlgError" in str(e):
                raise SingularMatrix(str(e))
            raise NoConvergenceError(str(e))
        
        x = x_res
        
        # Insert reference node voltage
        return self.toolkit.concatenate((x[:self.irefnode], self.toolkit.array([0.0]), x[self.irefnode:]))
    

    
    def get_diff(self, q, C, method=None):
        """Method used to calculate time derivative for charge storing elements (i_eq and g_eq)."""
        # Determine the active integrator based on step size variations
        h_last = getattr(self, '_dt_last', self._dt)
        self.active_integrator = self.base_integrator.check_order_drop(
            self._dt, h_last, getattr(self, '_is_first_step', False)
        )
        
        iq, geq = self.active_integrator.compute_derivatives(
            q_curr=q, C_curr=C, h_curr=self._dt, 
            q_last=self._qlast, iq_last=self._iqlast, h_last=h_last,
            is_first_step=getattr(self, '_is_first_step', False),
            toolkit=self.toolkit
        )
        
        self._iq = iq 
        self._effective_method = 'euler' if self.active_integrator.__class__.__name__ == 'EulerIntegrator' else self.par.method
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
        
        ## Respect a step controller injected by the caller (e.g. PIController);
        ## only fall back to the default IntegralController when none was set.
        if getattr(self, 'step_controller', None) is None:
            from pycircuit.circuit.stepcontroller import IntegralController
            self.step_controller = IntegralController()

        X = []
        self.irefnode=self.cir.get_node_index(refnode)
        n = self.cir.n
        if x0 is None:
            if self.par.uic:
                # Use Initial Conditions = skip DC operating point, start at zero
                x0 = self.toolkit.zeros(n)
            else:
                # Use DC operating point if x0 is not provided
                from pycircuit.circuit.dcanalysis import DC
                dc = DC(self.cir)
                try:
                    dc_res = dc.solve()
                    x0 = dc_res.x
                except Exception:
                    x0 = self.toolkit.zeros(n)
            x = x0
        else:
            x = x0 
        
        self.base_integrator = self._get_integrator()
        hist_len = max(2, self.base_integrator.get_required_history())
        q0 = self.cir.q(x)
        self._qlast = self.toolkit.array([q0 for _ in range(hist_len)])
        self._iqlast = self.toolkit.zeros((hist_len, n))
        
        X.append(copy(x))
        if hasattr(self.cir, 'accept_step'):
            self.cir.accept_step(0.0, X[-1], self.epar)
        
        timelist = []
        self._is_first_step = True
        t = 0.0
        max_step = timestep
        dt = timestep
        TRTOL = 7.0
        ## Bound the number of consecutive LTE rejections at a single time point.
        ## Near a source discontinuity (e.g. a VPulse corner) the truncation-error
        ## estimate can stay above tolerance for arbitrarily small steps while the
        ## stored history is frozen; without a cap the step size collapses and the
        ## solve grinds indefinitely.  After MAX_REJECT rejections we accept the
        ## already-converged Newton solution (only its LTE was too high) so time
        ## advances and the integrator history refreshes.
        reject_count = 0
        MAX_REJECT = 3
        
        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))
        abstol = self.toolkit.concatenate((self.par.iabstol * ones_nodes,
                                         self.par.vabstol * ones_branches))

        was_break_step = False
        while t < tend:
            # --- BREAKPOINT HANDLING ---
            # A breakpoint generally signifies a mathematical discontinuity (such as 
            # the sharp corner of a VPulse/IPulse square wave).
            # 
            # When the solver approaches a breakpoint, it does three vital things to 
            # maintain mathematical stability:
            # 1. Truncates Step (dt): If normal step size overshoots the breakpoint, 
            #    it forces dt to land *exactly* on the breakpoint timestamp.
            # 2. Flags the Breakpoint: was_break_step is set to True.
            # 3. Resets Integrator History: Immediately after crossing the breakpoint 
            #    (was_break_step == True in the next iteration), it sets 
            #    `self._is_first_step = True`. 
            #    
            # Why reset the history? Integrators (like Gear2 or Trapezoidal) use 
            # past state history to fit a smooth mathematical polynomial. If they 
            # tried to fit a polynomial across a sharp discontinuous edge, the 
            # simulation would suffer from massive artificial ringing and overshoot.
            # Resetting the history forces a drop to a safer 1st-order method 
            # (like Backward Euler) to gracefully navigate the corner and rebuild.
            if was_break_step:
                self._is_first_step = True
            
            next_t_break = self.cir.next_event(t)
            
            # Ensure next_t_break strictly advances time to avoid infinite dt=0 loops
            if next_t_break <= t + self.par.minbreak * max(abs(t), 1.0):
                next_t_break = self.cir.next_event(t + (self.par.minbreak * 1e3) * max(abs(t), 1.0))
            
            if t + dt > next_t_break:
                dt = float(next_t_break - t)
                was_break_step = True
            else:
                was_break_step = False
                
            if t + dt > tend:
                dt = tend - t
            
            self._dt = dt
            next_t = t + dt
            
            try:
                x, feval, J, f = self.solve_timestep(X[-1], next_t, provided_function=provided_function)
            except NoConvergenceError:
                dt = dt * 0.25
                if dt < getattr(self.par, 'minstep', 1e-18):
                    raise RuntimeError(f"Transient solver failed to converge: timestep shrank below {getattr(self.par, 'minstep', 1e-18):g}s at t={t}")
                continue
                
            if not fixed_timestep:
                accept, dt_next = self.step_controller.evaluate_step(
                    x_curr=x,
                    x_last=X[-1],
                    q_curr=self.cir.q(x),
                    q_last_hist=self._qlast,
                    iq_last_hist=self._iqlast,
                    h_curr=dt,
                    h_last=getattr(self, '_dt_last', dt),
                    is_first_step=self._is_first_step,
                    J=J,
                    active_integrator=self.active_integrator,
                    irefnode=self.irefnode,
                    reltol=self.par.reltol,
                    abstol=abstol,
                    toolkit=self.toolkit,
                    max_step=max_step,
                    TRTOL=TRTOL
                )
                
                if not accept and reject_count < MAX_REJECT:
                    reject_count += 1
                    dt = dt_next
                    if dt < getattr(self.par, 'minstep', 1e-18):
                        raise RuntimeError(f"Transient solver integration error: timestep shrank below {getattr(self.par, 'minstep', 1e-18):g}s at t={t}")
                    continue
                elif not accept:
                    ## Rejection cap reached: accept the converged solution at the
                    ## current (small) step, but grow dt back toward max_step so we
                    ## escape the collapsed regime instead of crawling step by step.
                    next_dt = min(max_step, dt * 10.0)
                else:
                    next_dt = dt_next
                reject_count = 0

            t = next_t
            timelist.append(t)
            X.append(copy(x))
            
            if hasattr(self.cir, 'accept_step'):
                self.cir.accept_step(t, X[-1], self.epar)
            
            # --- INTEGRATOR HISTORY RING BUFFERS ---
            # To support 2nd-order (and higher) integration methods, we must preserve the 
            # charge (q) and current/derivative (iq) of previous timesteps.
            # We push the newest values to index 0, and slice off the oldest `[:-1]` to 
            # maintain a constant buffer size (e.g. size 2 for Gear2).
            # This acts as a mathematical sliding window across the simulation time.
            self._iqlast = self.toolkit.concatenate((self.toolkit.array([self._iq]), self._iqlast))[:-1]
            self._qlast = self.toolkit.concatenate((self.toolkit.array([self.cir.q(x)]), self._qlast))[:-1]
            self._dt_last = dt
            
            self._is_first_step = False
            
            if not fixed_timestep:
                dt = next_dt
            
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
        self.base_integrator = self._get_integrator()
        hist_len = max(2, self.base_integrator.get_required_history())
        q0 = self.cir.q(x)
        self._qlast = self.toolkit.array([q0 for _ in range(hist_len)])
        self._iqlast = self.toolkit.zeros((hist_len, n))
        
        X.append(copy(x))
        if hasattr(self.cir, 'accept_step'):
            self.cir.accept_step(0.0, X[-1], self.epar)
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
            
            # Prepare eval_FJ callback for SchurCoupledNewton
            def eval_FJ(xr, h_curr):
                x_full = self.toolkit.zeros(n)
                x_full[:self.irefnode] = xr[:self.irefnode]
                x_full[self.irefnode+1:] = xr[self.irefnode:]
                
                self._dt = h_curr
                C = self.cir.C(x_full)
                q = self.cir.q(x_full)
                iq, Geq = self.get_diff(q, C)
                F = self.cir.i(x_full) + iq + self.cir.u(t + h_curr, analysis=self.par.analysis)
                J = self.cir.G(x_full) + Geq
                
                F_r = self.toolkit.delete(F, self.irefnode)
                J_r = self.toolkit.delete(J, self.irefnode, axis=0)
                J_r = self.toolkit.delete(J_r, self.irefnode, axis=1)
                
                eps = max(1e-8 * h_curr, 1e-15)
                self._dt = h_curr + eps
                iq_p, _ = self.get_diff(q, C)
                F_p = self.cir.i(x_full) + iq_p + self.cir.u(t + h_curr + eps, analysis=self.par.analysis)
                
                self._dt = h_curr - eps
                iq_m, _ = self.get_diff(q, C)
                F_m = self.cir.i(x_full) + iq_m + self.cir.u(t + h_curr - eps, analysis=self.par.analysis)
                
                J_h = (F_p - F_m) / (2 * eps)
                J_h_r = self.toolkit.delete(J_h, self.irefnode)
                self._dt = h_curr
                
                def calc_E(x_val, h_val):
                    if self._is_first_step: return 0.0
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
                
                E = calc_E(x_full, h_curr)
                
                gamma_min, gamma_max = 0.7, 3.0
                if analytical_eh and ((gamma_min - 1.0) * TRTOL <= E <= (gamma_max - 1.0) * TRTOL) and not self._is_first_step:
                    E_x_r = self.toolkit.zeros(len(F_r))
                    E_h = 1.0
                    E = 0.0
                elif self._is_first_step:
                    E_x_r = self.toolkit.zeros(len(F_r))
                    E_h = 1.0
                else:
                    if analytical_eh:
                        p = 3.0 if self._effective_method in ("trapezoidal", "trap") else 2.0
                        E_h = p * (E + TRTOL) / h_curr
                    else:
                        E_h = (calc_E(x_full, h_curr + eps) - calc_E(x_full, h_curr - eps)) / (2 * eps)
                    
                    E_x_r = self.toolkit.zeros(len(F_r))
                    if not analytical_eh:
                        for i in range(len(F_r)):
                            idx = i if i < self.irefnode else i + 1
                            x_p = copy(x_full)
                            x_p[idx] += eps
                            x_m = copy(x_full)
                            x_m[idx] -= eps
                            E_x_r[i] = (calc_E(x_p, h_curr) - calc_E(x_m, h_curr)) / (2 * eps)
                            
                return F_r, J_r, J_h_r, E, E_x_r, E_h

            def limiter_func(xr_next, xr_curr):
                x_next_full = self.toolkit.insert(xr_next, self.irefnode, 0.0)
                x_curr_full = self.toolkit.insert(xr_curr, self.irefnode, 0.0)
                x_next_full = self.cir.limit(x_next_full, x_curr_full, self.epar)
                return self.toolkit.concatenate((x_next_full[:self.irefnode], x_next_full[self.irefnode+1:]))

            from pycircuit.circuit.nrsolver import SchurCoupledNewton, NoConvergenceError
            solver = SchurCoupledNewton()
            
            x_curr_r = self.toolkit.concatenate((X[-1][:self.irefnode], X[-1][self.irefnode+1:]))
            abstol_r = self.toolkit.concatenate((abstol[:self.irefnode], abstol[self.irefnode+1:]))
            
            try:
                (x_next_r, h_next), _ = solver.solve_system(
                    (x_curr_r, h),
                    eval_FJ,
                    self.toolkit,
                    reltol,
                    abstol_r,
                    xtol,
                    self.par.maxiter,
                    limiter=limiter_func
                )
                converged = True
            except NoConvergenceError:
                converged = False
                
            if converged:
                x_next = self.toolkit.zeros(n)
                x_next[:self.irefnode] = x_next_r[:self.irefnode]
                x_next[self.irefnode+1:] = x_next_r[self.irefnode:]
                
                h = h_next
                x_curr = x_next
                    
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
            
            if hasattr(self.cir, 'accept_step'):
                self.cir.accept_step(t, X[-1], self.epar)
            
            self._dt = h_curr
            self._dt_last = h_curr
            self._is_first_step = False
            self._iqlast = self.toolkit.concatenate((self.toolkit.array([self._iq]), self._iqlast))[:-1]
            self._qlast = self.toolkit.concatenate((self.toolkit.array([self.cir.q(x_curr)]), self._qlast))[:-1]
            
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
