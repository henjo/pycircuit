import sys
sys.path.append("/home/andreas/sources/pycircuit")
import numpy as np
from copy import copy

import pycircuit.circuit._numeric as numeric
from pycircuit.circuit.transient import Transient
from pycircuit.circuit import SubCircuit, VS, R, C, gnd

c = SubCircuit()
c['VS'] = VS(1, gnd, v=10)
c['R1'] = R(1, 2, r=10)
c['C1'] = C(2, gnd, c=1e-9)
tran = Transient(c, toolkit=numeric)

# Test function
def _solve_coupled(tran, tend=1e-6, timestep=1e-8):
    tran.irefnode = tran.cir.get_node_index(gnd)
    n = tran.cir.n
    x0 = tran.toolkit.zeros(n)
    
    a,b,b_ = tran._method[tran.par.method]
    q0 = tran.cir.q(x0)
    tran._qlast = tran.toolkit.array([q0 for _ in range(max(2, len(a)))])
    tran._iqlast = tran.toolkit.zeros((max(2, len(b)), n))
    
    X = [copy(x0)]
    timelist = []
    
    tran._is_first_step = True
    t = 0.0
    h = timestep
    TRTOL = 7.0
    
    ones_nodes = tran.toolkit.ones(len(tran.cir.nodes))
    ones_branches = tran.toolkit.ones(len(tran.cir.branches))
    abstol = tran.toolkit.concatenate((tran.par.iabstol * ones_nodes,
                                     tran.par.vabstol * ones_branches))
    reltol = tran.par.reltol
    xtol = 1e-6
    
    while t < tend:
        if t + h > tend:
            h = tend - t
            
        x_curr = copy(X[-1])
        h_curr = h
        converged = False
        
        for iter_idx in range(tran.par.maxiter):
            # Compute F and J
            tran._dt = h_curr
            C = tran.cir.C(x_curr)
            q = tran.cir.q(x_curr)
            iq, Geq = tran.get_diff(q, C)
            F = tran.cir.i(x_curr) + iq + tran.cir.u(t + h_curr, analysis=tran.par.analysis)
            J = tran.cir.G(x_curr) + Geq
            
            F_r = tran.toolkit.delete(F, tran.irefnode)
            J_r = tran.toolkit.delete(J, tran.irefnode, axis=0)
            J_r = tran.toolkit.delete(J_r, tran.irefnode, axis=1)
            
            # Finite difference J_h
            eps = max(1e-8 * h_curr, 1e-15)
            tran._dt = h_curr + eps
            iq_p, Geq_p = tran.get_diff(q, C)
            F_p = tran.cir.i(x_curr) + iq_p + tran.cir.u(t + h_curr + eps, analysis=tran.par.analysis)
            
            tran._dt = h_curr - eps
            iq_m, Geq_m = tran.get_diff(q, C)
            F_m = tran.cir.i(x_curr) + iq_m + tran.cir.u(t + h_curr - eps, analysis=tran.par.analysis)
            
            J_h = (F_p - F_m) / (2 * eps)
            J_h_r = tran.toolkit.delete(J_h, tran.irefnode)
            tran._dt = h_curr
            
            # Compute E, E_x, E_h
            def calc_E(x_val, h_val):
                if tran._is_first_step:
                    return 0.0
                q_val = tran.cir.q(x_val)
                if tran.par.method == "trapezoidal":
                    dd2 = (q_val - tran._qlast[0]) / h_val - tran._iqlast[0]
                    dd2 = dd2 * 2.0 / h_val
                    lte = 1.0 / 12.0 * (h_val**3) * dd2
                elif tran.par.method == "gear2":
                    dd1_n = (q_val - tran._qlast[0]) / h_val
                    dd1_nm1 = (tran._qlast[0] - tran._qlast[1]) / tran._dt_last
                    dd2_n = (dd1_n - dd1_nm1) / (h_val + tran._dt_last)
                    lte = (h_val**2) * (h_val + tran._dt_last) / 3.0 * dd2_n
                else:
                    lte = tran.toolkit.zeros(n)
                etol = reltol * tran.toolkit.maximum(tran.toolkit.abs(q_val), tran.toolkit.abs(tran._qlast[0])) + abstol
                return tran.toolkit.max(tran.toolkit.abs(lte) / etol) - TRTOL
            
            E = calc_E(x_curr, h_curr)
            
            if tran._is_first_step:
                E_x_r = tran.toolkit.zeros(len(F_r))
                E_h = 1.0 # arbitrary non-zero to avoid division by zero
            else:
                E_p = calc_E(x_curr, h_curr + eps)
                E_m = calc_E(x_curr, h_curr - eps)
                E_h = (E_p - E_m) / (2 * eps)
                
                E_x_r = tran.toolkit.zeros(len(F_r))
                for i in range(len(F_r)):
                    # map i to actual index
                    idx = i if i < tran.irefnode else i + 1
                    x_p = copy(x_curr)
                    x_p[idx] += eps
                    E_xp = calc_E(x_p, h_curr)
                    
                    x_m = copy(x_curr)
                    x_m[idx] -= eps
                    E_xm = calc_E(x_m, h_curr)
                    
                    E_x_r[i] = (E_xp - E_xm) / (2 * eps)
            
            # Schur complement solve
            rhs = np.column_stack([-F_r, -J_h_r])
            dx_res = tran.toolkit.linearsolver(J_r, rhs)
            
            dx_0 = dx_res[:, 0]
            dx_h = dx_res[:, 1]
            
            if tran._is_first_step:
                dh = 0.0
            else:
                denom = tran.toolkit.dot(E_x_r, dx_h) + E_h
                if abs(denom) < 1e-20:
                    dh = 0.0
                else:
                    dh = (-E - tran.toolkit.dot(E_x_r, dx_0)) / denom
            
            # limit dh to prevent negative or overly large steps
            dh = max(-0.9 * h_curr, min(2.0 * h_curr, dh))
            
            dx_r = dx_0 + dx_h * dh
            
            x_next = copy(x_curr)
            x_next[:tran.irefnode] += dx_r[:tran.irefnode]
            x_next[tran.irefnode+1:] += dx_r[tran.irefnode:]
            
            x_curr = x_next
            h_curr += dh
            
            if tran.toolkit.all(tran.toolkit.abs(dx_r) < reltol * tran.toolkit.maximum(tran.toolkit.abs(x_curr[x_curr != 0]), 1e-12) + xtol):
                converged = True
                break
                
        if not converged:
            print("Failed to converge at t =", t)
            break
            
        t += h_curr
        h = h_curr
        X.append(copy(x_curr))
        timelist.append(t)
        
        tran._dt = h_curr
        tran._dt_last = h_curr
        tran._is_first_step = False
        tran._iqlast = tran.toolkit.concatenate((tran.toolkit.array([iq]), tran._iqlast))[:-1]
        tran._qlast = tran.toolkit.concatenate((tran.toolkit.array([q]), tran._qlast))[:-1]
        
    print("Coupled simulation completed, steps:", len(timelist))

_solve_coupled(tran)
