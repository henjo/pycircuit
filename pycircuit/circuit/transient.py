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
        lte_formula = getattr(self.par, 'lte_formula', 'classic')
        if self.par.method == 'euler':
            return EulerIntegrator(lte_formula)
        elif self.par.method in ('trap', 'trapezoidal'):
            return TrapezoidalIntegrator(lte_formula)
        elif self.par.method == 'gear2':
            return Gear2Integrator(lte_formula)
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
         Parameter(name='lte_formula',
                   desc="Local truncation error formula: 'ywr' "
                        "(Yao-Wang-Roychowdhury DAE LTE, ICECS 2014) or 'classic'",
                   unit='',
                   default='ywr'),
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
        max_step = timestep
        TRTOL = 7.0
        minstep = getattr(self.par, 'minstep', 1e-18)

        ones_nodes = self.toolkit.ones(len(self.cir.nodes))
        ones_branches = self.toolkit.ones(len(self.cir.branches))
        abstol = self.toolkit.concatenate((self.par.iabstol * ones_nodes,
                                         self.par.vabstol * ones_branches))
        reltol = self.par.reltol

        ## Coupled adaptive time-stepping, Fang, "A New Time-Stepping Method for
        ## Circuit Simulation" (DAC 2013).  The circuit solution and the step size
        ## are co-determined at each time point: converge the circuit at a fixed
        ## step (stage 1), then bring the local truncation error into the accept
        ## band by driving the step with Gear's formula and re-solving.  This is
        ## the paper's robust "approximate Newton" (sec. 3.4); it replaces the
        ## exact (N+1) Schur update, which is very sensitive to step changes and
        ## collapses the step size.  The LTE evaluation and Gear prediction are
        ## shared with the standard adaptive controller (IntegralController) so
        ## the coupled and adaptive paths stay consistent.
        from pycircuit.circuit.stepcontroller import IntegralController
        from pycircuit.circuit.nrsolver import NoConvergenceError
        controller = IntegralController()
        MAX_LTE_ITERS = 10

        while t < tend:
            if t + h > tend:
                h = tend - t

            ## Co-determine (x, h): converge the circuit at h_curr, evaluate the
            ## LTE, and while it is above the band shrink the step (Gear-predicted
            ## by the controller) and re-solve.
            h_curr = h
            x_curr = copy(X[-1])
            h_next = h
            for lte_iter in range(MAX_LTE_ITERS):
                self._dt = h_curr
                try:
                    x_new, feval, J, f = self.solve_timestep(
                        X[-1], t + h_curr, provided_function=provided_function)
                except NoConvergenceError:
                    ## Circuit did not converge at this step -> shrink and retry.
                    h_curr *= 0.25
                    if h_curr < minstep:
                        raise RuntimeError(f"Coupled transient: timestep shrank below {minstep:g}s at t={t}")
                    continue

                accept, h_next = controller.evaluate_step(
                    x_curr=x_new, x_last=X[-1],
                    q_curr=self.cir.q(x_new),
                    q_last_hist=self._qlast, iq_last_hist=self._iqlast,
                    h_curr=h_curr, h_last=getattr(self, '_dt_last', h_curr),
                    is_first_step=self._is_first_step, J=J,
                    active_integrator=self.active_integrator,
                    irefnode=self.irefnode, reltol=reltol, abstol=abstol,
                    toolkit=self.toolkit, max_step=max_step, TRTOL=TRTOL)

                x_curr = x_new
                if accept or lte_iter == MAX_LTE_ITERS - 1:
                    break
                if h_next < minstep:
                    raise RuntimeError(f"Coupled transient: timestep shrank below {minstep:g}s at t={t}")
                h_curr = h_next

            t += h_curr
            timelist.append(t)
            X.append(copy(x_curr))

            if hasattr(self.cir, 'accept_step'):
                self.cir.accept_step(t, X[-1], self.epar)

            self._dt = h_curr
            self._dt_last = h_curr
            self._is_first_step = False
            self._iqlast = self.toolkit.concatenate((self.toolkit.array([self._iq]), self._iqlast))[:-1]
            self._qlast = self.toolkit.concatenate((self.toolkit.array([self.cir.q(x_curr)]), self._qlast))[:-1]

            ## Next step: Gear-predicted size (already bounded by max_step).
            h = min(max_step, max(h_next, minstep))
            
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
