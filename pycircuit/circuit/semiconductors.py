# -*- coding: utf-8 -*-
from pycircuit.circuit import Circuit, Parameter, defaultepar

class Semiconductor(Circuit):
    """Base class for non-linear semiconductors with automatic Jacobians"""
    
    def G(self, x, epar=defaultepar):
        try:
            import jax
            if isinstance(x, jax.numpy.ndarray) or hasattr(self.toolkit, 'add_at'):
                # We can use JAX
                jac_func = jax.jacfwd(lambda v: self.i(v, epar))
                return self.toolkit.array(jac_func(x))
        except ImportError:
            pass
            
        # Fallback to central difference
        n = self.n
        G_mat = self.toolkit.zeros((n, n))
        eps = 1e-6
        for j in range(n):
            import copy
            x_plus = list(x)
            x_minus = list(x)
            x_plus[j] += eps
            x_minus[j] -= eps
            
            i_plus = self.i(self.toolkit.array(x_plus), epar)
            i_minus = self.i(self.toolkit.array(x_minus), epar)
            
            for k in range(n):
                G_mat[k, j] = (i_plus[k] - i_minus[k]) / (2 * eps)
        return G_mat

    def C(self, x, epar=defaultepar):
        # Default C is zero unless q is overridden
        if not hasattr(self.__class__, 'q') or self.__class__.q is Circuit.q:
            return self.toolkit.zeros((self.n, self.n))
            
        try:
            import jax
            if isinstance(x, jax.numpy.ndarray) or hasattr(self.toolkit, 'add_at'):
                jac_func = jax.jacfwd(lambda v: self.q(v, epar))
                return self.toolkit.array(jac_func(x))
        except ImportError:
            pass
            
        # Fallback to central difference
        n = self.n
        C_mat = self.toolkit.zeros((n, n))
        eps = 1e-6
        for j in range(n):
            x_plus = list(x)
            x_minus = list(x)
            x_plus[j] += eps
            x_minus[j] -= eps
            
            q_plus = self.q(self.toolkit.array(x_plus), epar)
            q_minus = self.q(self.toolkit.array(x_minus), epar)
            
            for k in range(n):
                C_mat[k, j] = (q_plus[k] - q_minus[k]) / (2 * eps)
        return C_mat

class BJT(Semiconductor):
    """NPN Bipolar Junction Transistor (Ebers-Moll Level 1)"""
    terminals = ('c', 'b', 'e')
    instparams = [
        Parameter('IS', 'Saturation current', default=1e-14, unit='A'),
        Parameter('BF', 'Forward beta', default=100.0),
        Parameter('BR', 'Reverse beta', default=1.0),
        Parameter('VA', 'Forward Early voltage', default=100.0, unit='V'),
    ]

    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        v_c = x[0]
        v_b = x[1]
        v_e = x[2]
        
        v_be = v_b - v_e
        v_bc = v_b - v_c
        v_ce = v_c - v_e
        
        VT = toolkit.kboltzmann * epar.T / toolkit.qelectron
        IS = params.get('IS', 1e-14)
        BF = params.get('BF', 100.0)
        BR = params.get('BR', 1.0)
        VA = params.get('VA', 100.0)
        
        # --- EBERS-MOLL (LEVEL 1) BJT MODEL ---
        # The BJT is modeled as two coupled diodes (Base-Emitter and Base-Collector).
        # These diodes share the base region and create transistor action.
        i_be_diode = IS * (toolkit.exp(v_be / VT) - 1.0)
        i_bc_diode = IS * (toolkit.exp(v_bc / VT) - 1.0)
        
        # Early effect: models the modulation of the effective base width by the
        # collector-emitter voltage, which causes output conductance to increase.
        early_factor = (1.0 + v_ce / VA) if VA > 0 else 1.0
        
        # The transport current (i_cc) is the minority carrier current injected 
        # completely across the base.
        i_cc = (i_be_diode - i_bc_diode) * early_factor
        
        # KCL: Base current consists of the forward and reverse recombination losses.
        i_b_val = (i_be_diode / BF) + (i_bc_diode / BR)
        
        # Collector current is the transport current minus reverse recombination.
        i_c_val = i_cc - (i_bc_diode / BR)
        
        # Emitter current is exactly what goes in from C and B (KCL check: Ic+Ib+Ie=0)
        i_e_val = -i_cc - (i_be_diode / BF)
        
        return toolkit.array([i_c_val, i_b_val, i_e_val])

    def i(self, x, epar=defaultepar):
        params = {
            'IS': self.iparv.IS,
            'BF': self.iparv.BF,
            'BR': self.iparv.BR,
            'VA': self.iparv.VA
        }
        return self.eval_i_pure(x, params, epar, self.toolkit)


class JFET(Semiconductor):
    """N-Channel JFET (Shichman-Hodges)"""
    terminals = ('d', 'g', 's')
    instparams = [
        Parameter('VTO', 'Threshold voltage', default=-2.0, unit='V'),
        Parameter('BETA', 'Transconductance parameter', default=1e-4, unit='A/V^2'),
        Parameter('LAMBDA', 'Channel length modulation', default=0.0, unit='1/V'),
        Parameter('IS', 'Gate p-n saturation current', default=1e-14, unit='A'),
    ]

    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        v_d = x[0]
        v_g = x[1]
        v_s = x[2]
        
        v_gs = v_g - v_s
        v_ds = v_d - v_s
        v_gd = v_g - v_d
        
        where = getattr(toolkit, 'where', None)
        if where is None:
            import jax.numpy as jnp
            where = jnp.where
            
        VTO = params.get('VTO', -2.0)
        BETA = params.get('BETA', 1e-4)
        LAMBDA = params.get('LAMBDA', 0.0)
        IS = params.get('IS', 1e-14)
        
        # --- SHICHMAN-HODGES JFET MATHEMATICAL REGIONS ---
        # The JFET has 3 distinct operating regions (Cutoff, Triode/Linear, Saturation).
        # To vectorize this effectively on a GPU without standard Python `if` statements, 
        # we use branchless mathematical masking (`where`).

        # Symmetric Device Swap: If V_ds < 0, the drain and source swap roles.
        is_reverse = v_ds < 0
        v_ds_eff = where(is_reverse, -v_ds, v_ds)
        v_gs_eff = where(is_reverse, v_gd, v_gs)
        
        # Region definitions:
        # 1. Cutoff: V_gs is more negative than the threshold (VTO)
        # 2. Saturation: V_ds is greater than the overdrive voltage (V_gs - VTO)
        # 3. Triode: V_ds is smaller than the overdrive voltage
        v_ds_sat = v_gs_eff - VTO
        is_cutoff = v_gs_eff <= VTO
        is_sat = v_ds_eff >= v_ds_sat
        
        # Triode formula: I_d = Beta * [2*(V_gs - VTO)*V_ds - V_ds^2] * (1 + Lambda*V_ds)
        i_d_triode = BETA * v_ds_eff * (2.0 * (v_gs_eff - VTO) - v_ds_eff) * (1.0 + LAMBDA * v_ds_eff)
        i_d_sat = BETA * (v_gs_eff - VTO)**2 * (1.0 + LAMBDA * v_ds_eff)
        
        i_d_fwd = where(is_cutoff, 0.0, where(is_sat, i_d_sat, i_d_triode))
        i_d = where(is_reverse, -i_d_fwd, i_d_fwd)
        
        VT = toolkit.kboltzmann * epar.T / toolkit.qelectron
        i_gs = IS * (toolkit.exp(v_gs / VT) - 1.0)
        i_gd = IS * (toolkit.exp(v_gd / VT) - 1.0)
        
        i_g = i_gs + i_gd
        i_drain = i_d - i_gd
        i_source = -i_d - i_gs
        
        return toolkit.array([i_drain, i_g, i_source])

    def i(self, x, epar=defaultepar):
        params = {
            'VTO': self.iparv.VTO,
            'BETA': self.iparv.BETA,
            'LAMBDA': self.iparv.LAMBDA,
            'IS': self.iparv.IS
        }
        return self.eval_i_pure(x, params, epar, self.toolkit)

class ZenerDiode(Semiconductor):
    """Diode with Avalanche Breakdown"""
    terminals = ('plus', 'minus')
    instparams = [
        Parameter('IS', 'Saturation current', default=1e-14, unit='A'),
        Parameter('BV', 'Reverse breakdown voltage', default=5.0, unit='V'),
        Parameter('IBV', 'Reverse breakdown current', default=1e-10, unit='A'),
    ]

    def i(self, x, epar=defaultepar):
        v = x[0] - x[1]
        VT = self.toolkit.kboltzmann * epar.T / self.toolkit.qelectron
        
        i_fwd = self.iparv.IS * (self.toolkit.exp(v / VT) - 1.0)
        
        v_zener = -v - self.iparv.BV
        i_zener = -self.iparv.IBV * (self.toolkit.exp(v_zener / VT) - 1.0)
        
        i_tot = i_fwd + i_zener
        return self.toolkit.array([i_tot, -i_tot])

class Varactor(Semiconductor):
    """Voltage-dependent Capacitor"""
    terminals = ('plus', 'minus')
    instparams = [
        Parameter('CJ0', 'Zero-bias junction capacitance', default=1e-12, unit='F'),
        Parameter('VJ', 'Junction potential', default=1.0, unit='V'),
        Parameter('M', 'Grading coefficient', default=0.5),
    ]

    def i(self, x, epar=defaultepar):
        return self.toolkit.zeros(2)
        
    def q(self, x, epar=defaultepar):
        v = x[0] - x[1]
        
        CJ0 = self.iparv.CJ0
        VJ = self.iparv.VJ
        M = self.iparv.M
        
        # JAX-friendly min function
        try:
            import jax.numpy as jnp
            is_tracer = isinstance(v, jnp.ndarray) or hasattr(v, 'aval')
        except ImportError:
            is_tracer = False
            
        if is_tracer or hasattr(self.toolkit, 'minimum'):
            minimum = getattr(self.toolkit, 'minimum', None)
            if minimum is None:
                import jax.numpy as jnp
                minimum = jnp.minimum
            v_eff = minimum(v, 0.99 * VJ)
        else:
            v_eff = v if v <= 0.99 * VJ else 0.99 * VJ
                
        q_val = CJ0 * VJ / (1.0 - M) * (1.0 - (1.0 - v_eff / VJ)**(1.0 - M))
        return self.toolkit.array([q_val, -q_val])
