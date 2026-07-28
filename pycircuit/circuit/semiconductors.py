# -*- coding: utf-8 -*-
from pycircuit.circuit import Circuit, Parameter, defaultepar

class Semiconductor(Circuit):
    """Base class for non-linear semiconductors with automatic Jacobians.

    These devices deliberately write no ``G``/``C`` of their own: the model
    lives once, in ``eval_i_pure``/``eval_q_pure``, and the Jacobian is
    obtained by differentiating it.  That is the invariant ``P9``/``P10`` in
    ``doc/architecture.md`` exist to protect -- a hand-written stamp that
    silently disagrees with ``i()`` still lets Newton converge, just to a
    slightly wrong answer.

    Which differentiation is used is the *toolkit's* decision, not this
    class's.  It used to be decided here by ``try: import jax``, i.e. by
    whether JAX happened to be installed rather than by which backend was
    running -- the same defect ``Toolkit.jacobian`` and ``Toolkit.derivative``
    were introduced to remove from ``elements.py``, which this module was
    never migrated onto.  With a symbolic toolkit that fell through to a
    central difference, perturbing a *sympy symbol* by ``1e-6``.
    """

    def _model_params(self):
        """Instance parameters as a plain dict, for the pure model functions."""
        return {p.name: getattr(self.iparv, p.name) for p in self.instparams}

    def G(self, x, epar=defaultepar):
        return self.toolkit.jacobian(self.eval_i_pure, x,
                                     self._model_params(), epar)

    def C(self, x, epar=defaultepar):
        ## No charge model means no capacitance -- nothing to differentiate.
        if not hasattr(self.__class__, 'q') or self.__class__.q is Circuit.q:
            return self.toolkit.zeros((self.n, self.n))

        return self.toolkit.jacobian(self.eval_q_pure, x,
                                     self._model_params(), epar)

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
        return self.eval_i_pure(x, self._model_params(), epar, self.toolkit)


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
        
        where = toolkit.where

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
        return self.eval_i_pure(x, self._model_params(), epar, self.toolkit)

class ZenerDiode(Semiconductor):
    """Diode with Avalanche Breakdown"""
    terminals = ('plus', 'minus')
    instparams = [
        Parameter('IS', 'Saturation current', default=1e-14, unit='A'),
        Parameter('BV', 'Reverse breakdown voltage', default=5.0, unit='V'),
        Parameter('IBV', 'Reverse breakdown current', default=1e-10, unit='A'),
    ]

    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        v = x[0] - x[1]
        VT = toolkit.kboltzmann * epar.T / toolkit.qelectron
        IS = params.get('IS', 1e-14)
        BV = params.get('BV', 5.0)
        IBV = params.get('IBV', 1e-10)
        
        i_fwd = IS * (toolkit.exp(v / VT) - 1.0)
        v_zener = -v - BV
        i_zener = -IBV * (toolkit.exp(v_zener / VT) - 1.0)
        
        i_tot = i_fwd + i_zener
        return toolkit.array([i_tot, -i_tot])

    def i(self, x, epar=defaultepar):
        return self.eval_i_pure(x, self._model_params(), epar, self.toolkit)

class Varactor(Semiconductor):
    """Voltage-dependent Capacitor"""
    terminals = ('plus', 'minus')
    instparams = [
        Parameter('CJ0', 'Zero-bias junction capacitance', default=1e-12, unit='F'),
        Parameter('VJ', 'Junction potential', default=1.0, unit='V'),
        Parameter('M', 'Grading coefficient', default=0.5),
    ]

    @staticmethod
    def eval_i_pure(x, params, epar, toolkit):
        return toolkit.zeros(2)

    @staticmethod
    def eval_q_pure(x, params, epar, toolkit):
        v = x[0] - x[1]
        CJ0 = params.get('CJ0', 1e-12)
        VJ = params.get('VJ', 1.0)
        M = params.get('M', 0.5)
        
        minimum = toolkit.minimum

        v_eff = minimum(v, 0.99 * VJ)
        q_val = CJ0 * VJ / (1.0 - M) * (1.0 - (1.0 - v_eff / VJ)**(1.0 - M))
        return toolkit.array([q_val, -q_val])

    def i(self, x, epar=defaultepar):
        return self.eval_i_pure(x, self._model_params(), epar, self.toolkit)

    def q(self, x, epar=defaultepar):
        return self.eval_q_pure(x, self._model_params(), epar, self.toolkit)
