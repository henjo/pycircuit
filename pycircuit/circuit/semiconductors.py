# -*- coding: utf-8 -*-
import numpy as np

from pycircuit.circuit import Circuit, Parameter, defaultepar

## STAGE 5(b) -- SPICE's EXPMAX treatment, and why it belongs in this module.
##
## The Ebers-Moll and Shockley expressions below all evaluate `exp(v/VT)`, and a
## Newton step can ask for a `v` these devices will never see in operation.  At
## T=300K, VT is 25.85 mV, so `exp()` overflows to `inf` above about 18.3 V --
## while the supply in the failing common-emitter case is 5 V and `exp(5/VT)` is
## perfectly finite.  **The excursion comes from the iteration, not the bias.**
##
## What overflow costs here is worse than a large number.  `inf - inf` in the
## transport current gives `nan`; `nan` reaches `x`; every later evaluation is
## `nan`.  Measured on that circuit: 12 of 405 model evaluations returned
## non-finite `i()` AND non-finite `G()`.
##
## It belongs on this base class rather than in each device because these models
## are differentiated by `toolkit.jacobian`.  Under a toolkit without autodiff
## that is a CENTRAL DIFFERENCE, so an overflow in `eval_i_pure` does not merely
## produce a large entry in `G` -- it produces `nan`, from `(inf - inf)/2h`.  One
## clamp at the source covers both the residual and the Jacobian.
##
## The linear extrapolation above the knee is deliberate, and is what SPICE does:
## clamping `exp` to a constant would zero the derivative and leave Newton with no
## gradient to descend, which trades a `nan` for a stall.  Continuing linearly
## keeps the function monotone and its derivative positive and finite.
EXP_ARG_MAX = 80.0


def _expl(toolkit, arg):
    """`exp(arg)`, extrapolated linearly above `EXP_ARG_MAX` instead of overflowing.

    Above the knee, `exp(a) ~ exp(m)(1 + (a - m))` -- value and first derivative
    both continuous at `a = m`, so Newton sees no discontinuity and still has a
    finite positive slope to work with.
    """
    lim = EXP_ARG_MAX
    e_lim = toolkit.exp(lim)
    try:
        if arg > lim:
            return e_lim * (1.0 + (arg - lim))
        return toolkit.exp(arg)
    except (TypeError, ValueError):
        ## Symbolic or array-valued argument: no comparison is possible, so hand
        ## back the plain exponential rather than guessing.  The clamp exists for
        ## the numeric Newton path, which is where the overflow happens.
        return toolkit.exp(arg)


def _depletion_charge(toolkit, v, CJ, VJ, M, FC=0.5):
    """SPICE's junction depletion charge, linearised above the `FC*VJ` knee.

    Below the knee this is the textbook expression; above it, SPICE continues
    with a quadratic whose value AND first derivative match at `v = FC*VJ`, so
    `C = dq/dv` stays finite, positive and increasing.

    **The alternative in this file is what stage 5+.3 exists to fix.**  `Varactor`
    freezes the charge with `min(v, 0.99*VJ)`, which makes `C = dq/dv` fall to
    EXACTLY ZERO in forward bias -- on the node with the largest physical
    capacitance in the circuit.  That is worse than inaccurate: it removes the
    state variable, and a Newton step that sees `C = 0` there takes a wildly wrong
    step.  One treatment, shared, so a new device cannot inherit the old defect.
    """
    knee = FC * VJ
    below = CJ * VJ / (1.0 - M) * (1.0 - (1.0 - v / VJ) ** (1.0 - M))
    try:
        if v <= knee:
            return below
    except (TypeError, ValueError):
        ## Symbolic or array-valued: no branch is possible, so use the smooth
        ## expression.  The knee only matters to the numeric Newton path.
        return below

    ## F1/F2/F3 are SPICE's constants; F1 is the charge AT the knee, so the two
    ## pieces meet there by construction rather than by choosing a scale factor.
    F1 = VJ / (1.0 - M) * (1.0 - (1.0 - FC) ** (1.0 - M))
    F2 = (1.0 - FC) ** (1.0 + M)
    F3 = 1.0 - FC * (1.0 + M)
    return CJ * (F1 + (1.0 / F2) * (F3 * (v - knee)
                                    + (M / (2.0 * VJ)) * (v * v - knee * knee)))

def _pnjlim(vnew, vold, VT, IS, toolkit):
    """SPICE's `pnjlim`: bound the per-iteration excursion of a junction voltage.

    Newton linearises an exponential, so a step taken from a lightly-biased point
    overshoots enormously -- and the model is then evaluated somewhere it has no
    business being.  `_expl` stops that becoming `nan`, but it does not stop the
    iteration wandering: with the clamp alone and no limiting, the common-emitter
    stage converges to a GENUINE but spurious operating point 200 V below the
    rails, with the base-collector junction forward biased at 0.91 V carrying the
    20 A the base resistor delivers.  That is a solution of the circuit equations;
    it is simply not the one anyone wants, and only limiting keeps Newton out of
    it.

    `vc` is the voltage at which the exponential's curvature starts to dominate.
    Above it the update is compressed logarithmically, which is what bounds the
    step without changing where the solution is -- the limiter only moves the
    point the next Jacobian is taken at, never the equations.
    """
    if IS <= 0.0:
        return vnew
    vc = VT * toolkit.log(VT / (IS * 1.414213562))
    if vnew > vc and vnew > 0.0:
        if vold > 0.0:
            arg = 1.0 + (vnew - vold) / VT
            return vold + VT * toolkit.log(arg) if arg > 0.0 else vc
        return VT * toolkit.log(vnew / VT)
    return vnew


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

    ## STAGE 5(a) -- WHICH JUNCTIONS THIS DEVICE HAS, AND WHICH TERMINAL MOVES.
    ##
    ## Each entry is `(anode, cathode, move)` as indices into `terminals`.  The
    ## third field is not redundant: a BJT's two junctions SHARE the base as anode,
    ## so limiting both by adjusting the anode would have the second undo the first.
    ## Moving the base for B-E and the collector for B-C keeps the two independent,
    ## which is the same decomposition the Ebers-Moll model itself uses.
    ##
    ## Empty by default: a device with no junction list gets no limiting, which is
    ## the right answer for anything whose i(v) is not an exponential.
    junctions = ()

    def limit(self, x, x0, epar=defaultepar):
        """Return a limited copy of `x` -- STATE-FREE, and that is the point.

        The convention here is that `limit` RETURNS the limited vector rather than
        mutating device state.  `Diode` does the opposite: it stores `_vlim` and
        linearises `G` around it, which works but hid a real bug for a long time --
        `SubCircuit.limit` passes a fancy-indexed COPY of `x` and discards the
        result, so a limiter that writes its argument has no effect at all, and the
        only limiter in the tree was the one that could not expose that.

        Returning the vector also makes limiting expressible on a traced backend,
        where device-private Python state cannot go.  Both `None` (the `Diode`
        convention) and a returned array are accepted by `SubCircuit.limit` during
        the transition.
        """
        if not self.junctions:
            return x

        VT = self.toolkit.kboltzmann * epar.T / self.toolkit.qelectron
        IS = getattr(self.iparv, 'IS', 0.0)
        try:
            out = np.array(x, dtype=float, copy=True)
        except (TypeError, ValueError):
            ## Symbolic x: limiting is a numeric Newton aid and has nothing to
            ## contribute here.  Returning it untouched is correct, not a fallback.
            return x

        x0a = np.asarray(x0, dtype=float)
        for anode, cathode, move in self.junctions:
            vnew = float(out[anode] - out[cathode])
            vold = float(x0a[anode] - x0a[cathode])
            vlim = _pnjlim(vnew, vold, VT, IS, self.toolkit)
            ## Reassign the moving terminal so the junction carries `vlim` exactly.
            if move == anode:
                out[anode] = out[cathode] + vlim
            else:
                out[cathode] = out[anode] - vlim
        return out

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
    ## B-E moves the base; B-C moves the collector.  See `junctions` on the base
    ## class for why the two must not both move the base.
    junctions = ((1, 2, 1), (1, 0, 0))
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
        ## `_expl`, not `exp` -- see EXP_ARG_MAX.  `inf - inf` here is what put
        ## `nan` into 12 of 405 evaluations on the failing common-emitter stage.
        i_be_diode = IS * (_expl(toolkit, v_be / VT) - 1.0)
        i_bc_diode = IS * (_expl(toolkit, v_bc / VT) - 1.0)
        
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

    ## STAGE 5+.1 -- THE CHARGE MODEL, and where these defaults come from.
    ##
    ## Review 0.1b ranked this above the missing limiter: `BJT` defined no `q`, so
    ## `Semiconductor.C` took its zero-matrix branch and the transistor had NO
    ## charge storage at all.  Measured before this: `C(x)` identically zero, and
    ## the same switching stage driven 1000x faster produced a **bit-identical**
    ## waveform -- there was no time constant in the device to respond with.  A
    ## bipolar transient was a sequence of DC solves, and every switching time it
    ## produced was wrong.
    ##
    ## SPICE defaults CJE/CJC/TF/TR to ZERO and expects a model card.  These do
    ## not, and the reason is consistency with this class rather than with SPICE:
    ## `BJT` already invents IS=1e-14, BF=100, VA=100 -- it is a usable default
    ## transistor, not a bare template.  Defaulting the charge to zero would leave
    ## the defect above fully in place for anyone who did not know to set CJE and
    ## TF, which is the confidently-wrong-answer class this work exists to remove.
    ##
    ## **These are a default model card, not a measurement.**  They are the round
    ## numbers of a generic small-signal NPN and are documented as such; any real
    ## comparison against silicon must supply its own.
    instparams = instparams + [
        Parameter('CJE', 'B-E zero-bias depletion capacitance', default=4.5e-12, unit='F'),
        Parameter('VJE', 'B-E junction potential', default=0.75, unit='V'),
        Parameter('MJE', 'B-E grading coefficient', default=0.33),
        Parameter('CJC', 'B-C zero-bias depletion capacitance', default=3.5e-12, unit='F'),
        Parameter('VJC', 'B-C junction potential', default=0.75, unit='V'),
        Parameter('MJC', 'B-C grading coefficient', default=0.33),
        Parameter('TF', 'Forward transit time', default=3e-10, unit='s'),
        Parameter('TR', 'Reverse transit time', default=4e-9, unit='s'),
        Parameter('FC', 'Forward-bias depletion coefficient', default=0.5),
    ]

    @staticmethod
    def eval_q_pure(x, params, epar, toolkit):
        """Terminal charges: depletion on both junctions plus diffusion via TF/TR.

        The two junctions both connect to the base, so the base carries the sum
        and each of the collector and emitter carries the negative of its own --
        which is what makes `sum(q) == 0`, i.e. the device stores no net charge.
        """
        v_c, v_b, v_e = x[0], x[1], x[2]
        v_be = v_b - v_e
        v_bc = v_b - v_c

        VT = toolkit.kboltzmann * epar.T / toolkit.qelectron
        IS = params.get('IS', 1e-14)
        FC = params.get('FC', 0.5)

        ## Depletion charge -- the same helper `Varactor` uses, so the FC
        ## linearisation exists once.
        q_dep_be = _depletion_charge(toolkit, v_be, params.get('CJE', 0.0),
                                     params.get('VJE', 0.75),
                                     params.get('MJE', 0.33), FC)
        q_dep_bc = _depletion_charge(toolkit, v_bc, params.get('CJC', 0.0),
                                     params.get('VJC', 0.75),
                                     params.get('MJC', 0.33), FC)

        ## Diffusion charge: minority carriers in transit, proportional to the
        ## junction current.  `_expl` for the same reason `eval_i_pure` uses it --
        ## this is differentiated by a central difference under a toolkit without
        ## autodiff, so an overflow here becomes `nan` in C, not a large entry.
        q_dif_be = params.get('TF', 0.0) * IS * (_expl(toolkit, v_be / VT) - 1.0)
        q_dif_bc = params.get('TR', 0.0) * IS * (_expl(toolkit, v_bc / VT) - 1.0)

        q_be = q_dep_be + q_dif_be
        q_bc = q_dep_bc + q_dif_bc

        return toolkit.array([-q_bc, q_be + q_bc, -q_be])

    def q(self, x, epar=defaultepar):
        return self.eval_q_pure(x, self._model_params(), epar, self.toolkit)


class JFET(Semiconductor):
    """N-Channel JFET (Shichman-Hodges)"""
    terminals = ('d', 'g', 's')
    ## terminals ('d', 'g', 's'): G-S moves the gate, G-D moves the drain.
    junctions = ((1, 2, 1), (1, 0, 0))
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
        i_gs = IS * (_expl(toolkit, v_gs / VT) - 1.0)
        i_gd = IS * (_expl(toolkit, v_gd / VT) - 1.0)
        
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
        
        i_fwd = IS * (_expl(toolkit, v / VT) - 1.0)
        v_zener = -v - BV
        ## The clamp is direction-agnostic, which is why the Zener gets it even
        ## though review 0.1b withheld `pnjlim` from it: this term is a second
        ## exponential running the OTHER way, and overflow there is just as fatal.
        i_zener = -IBV * (_expl(toolkit, v_zener / VT) - 1.0)
        
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
