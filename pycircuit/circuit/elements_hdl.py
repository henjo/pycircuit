# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Elements written in the Behavioural HDL -- the language, eating its own food.

Every element here is a re-expression of a hand-written ``elements.py``
element (or a well-known behavioral primitive) in contribution-statement
form, and ``tests/test_elements_hdl.py`` proves each one equivalent to its
hand-written counterpart -- stamp-for-stamp where the formulations match,
waveform-for-waveform where only the physics does (an HDL branch row may
legitimately be the hand-written row times -1; MNA does not care).

What this file demonstrates, capability by capability:

===============  ====================================================
element          HDL capability it exercises
===============  ====================================================
``RHdl``         I-contribution, static linear
``GHdl``         same, conductance form
``CHdl``         ``ddt`` -> charge routing
``LHdl``         V-contribution + self current probe + ``ddt`` of a
                 branch current (the inductor needs all three)
``VCCSHdl``      4 terminals, cross-branch probe
``VCVSHdl``      V-contribution controlled by a remote voltage
``DiodeHdl``     nonlinear static (``exp``), exact symbolic Jacobian
``VSinHdl``      time-dependent source terms via ``TIME`` (+ ``dudt``)
``IdtHdl``       the ``idt`` operator -> generated state equation
``IdtmodHdl``    ``idtmod`` -> state + floored wrap + generated
                 ``periodic_states`` (the Phase-2 gauge shift) + DC pin
``RThermalHdl``  self-heating: an extra node whose "voltage" is the
                 temperature rise, an RC to ambient, dissipated power
                 injected as a current, and the device's own resistance
                 reading the node back
``DiodeSpiceHdl``  the full SPICE level-1 junction -- series resistance
                 on an internal node with parameter-driven `Collapse`,
                 depletion + diffusion charge through ``ddt`` with the
                 ``fc`` linearisation, the ``is``/``vj``/``cj0``
                 temperature path, reverse breakdown, area scaling
``DiodeSpiceThermalHdl``  the same junction with the thermal node
                 bolted on -- four extra lines in ``analog()``
===============  ====================================================

Hand-written elements remain the defaults (``elements.py``): they carry
noise models, device limiting, and per-element tuning the HDL does not
generate yet.  This file is the proof of expressiveness and the natural
home for user macromodels; ``benchmarks/hdl_overhead.py`` measures what
the generated code costs against the hand-written stamps.
"""

import sympy

from pycircuit.utilities.param import Parameter
from pycircuit.circuit.hdl import (Behavioural, Branch, Node, Contribution,
                                   Collapse, ddt, idt, idtmod, TIME, TEMP)

## ---------------------------------------------------------------------
## ALIASED WITH A LEADING UNDERSCORE.  The reason it was NECESSARY is
## gone; the alias is kept because it is still the clearer thing to read.
##
## 2026-08-25: `generate_code` used to bind each class's parameter names
## by doing `analogfunc.__globals__.update(...)`, and `copy(cls.analog)`
## returns the *same* function object, so `__globals__` WAS this module's
## namespace.  `DiodeHdl` declares a parameter named `vt`, so importing
## this module replaced `elements_hdl.vt` with `Symbol('vt')` and any
## model defined below it that wrote `vt()` called a Symbol.
##
## That injection is now non-mutating (roadmap Phase 2): the parameter
## symbols go into a private copy of this module's globals, and nothing
## is rebound here.  `elements_hdl.vt` is no longer created at all, and
## `_vt is hdl.vt` -- both pinned by a test in test_hdl_diagnostics.py.
##
## The underscore is retained anyway: a helper shared by several classes
## reads better when the name it calls cannot be confused with any
## class's parameter, and `Parameter(name='_limexp')` is not something
## anyone writes.  It is now a readability choice, not a workaround.
## ---------------------------------------------------------------------
from pycircuit.circuit.hdl import vt as _vt
from pycircuit.circuit.hdl import var as _var
## `expl`, not `limexp`, and the reason is measured rather than
## stylistic.  `limexp` is deliberately NOT both-arms-safe: its discarded
## `exp(x)` overflows above 709 so that PCNR's shape detector can still
## read `exp(arg)` out of the expression and recover the junction's IS
## and VT.  That trade only pays if the model can actually participate
## in PCNR -- and `DiodeSpiceHdl` never can.
##
## 2026-08-25, CORRECTED: this comment used to say the disqualifier was
## the CHARGE, on the grounds that the layer refuses a participant that
## stores it.  That was true once and is not true now -- `pcnr.py:85`
## records the refusal being removed, precisely because it "cost every
## junction device with capacitance, which is to say every real one",
## and `test_pcnr_charge.py` proves the charge case.
##
## Measured instead, three elements differing by one thing each:
##
##     flat, no charge      chained=False   pcnr_i=True
##     flat, WITH charge    chained=False   pcnr_i=True
##     chained, with charge chained=True    pcnr_i=False
##
## The disqualifier is `var()`.  The shape detector reads `exp(arg)` out
## of the expression, and a let-chain hides it behind an intermediate
## symbol.  `DiodeSpiceHdl` is chained (and has a second branch for
## `rs`), so `hasattr(el, 'pcnr_i')` is False for every parameter set --
## the right answer, for a different reason than was written here.
##
## The consequence is larger than this one model: EVERY production model
## uses `var()`, so none of them can participate in PCNR.  That is the
## same root cause as the missing `pure_spec` in roadmap section 2 --
## the chained path costs two capabilities, not one.
##
## So the model paid the cost and got none of
## the benefit: 120 `overflow encountered in exp` warnings out of ONE
## 5 V DC solve of the self-heating variant, which makes it unusable
## under the `np.seterr(all='raise')` anyone debugging a convergence
## failure will set.
##
## The self-heating variant is where it bites, and the mechanism is
## worth recording: with a thermal node the exponential's own SCALE
## `n*k*T/q` is a solution variable, so `vd/nvt` is unbounded in a way
## the isothermal model's never is.  A Newton iterate that visits
## T ~ 1 K makes `0.7/nvt` about 7700.  Confirmed by raising the
## temperature floor to 200 K: 120 warnings become 0.
##
## `expl` is the same function below 230.26 -- exp, exactly -- and
## continues with exp's own third-order Taylor above it, with both arms
## clamped to their own side.  Measured across 5 V, 20 V and 100 V
## drives: identical solutions to the last digit, zero overflows.  It is
## also what JUNCAP200, this repo's own junction model, is written in.
from pycircuit.circuit.hdl import expl as _expl
from pycircuit.circuit.hdl import hypsmooth as _hypsmooth
from pycircuit.circuit.hdl import safe_abs as _safe_abs
from pycircuit.circuit.hdl import white_noise as _white_noise
from pycircuit.circuit.hdl import flicker_noise as _flicker_noise
from pycircuit.circuit.hdl import KBOLTZMANN as _KB
from pycircuit.circuit.hdl import QELECTRON as _QE


class RHdl(Behavioural):
    """Resistor: ``I(p,n) <+ V(p,n)/r``.  (No noise model -- use
    ``elements.R`` when noise matters.)"""
    instparams = [Parameter(name='r', desc='Resistance', unit='ohm',
                            default=1e3)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, b.V / r)


class GHdl(Behavioural):
    """Conductor: ``I(p,n) <+ g*V(p,n)``."""
    instparams = [Parameter(name='g', desc='Conductance', unit='S',
                            default=1e-3)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, g * b.V)


class CHdl(Behavioural):
    """Capacitor: ``I(p,n) <+ ddt(c*V(p,n))``."""
    instparams = [Parameter(name='c', desc='Capacitance', unit='F',
                            default=1e-12)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, ddt(c * b.V))


class LHdl(Behavioural):
    """Inductor: ``V(p,n) <+ ddt(L*I(p,n))``.

    The V-contribution creates the branch-current unknown; probing the
    same branch's current inside ``ddt`` closes the loop -- three language
    features one classical element needs at once.
    """
    instparams = [Parameter(name='L', desc='Inductance', unit='H',
                            default=1e-9)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.V, ddt(L * b.I))


class VCCSHdl(Behavioural):
    """Voltage-controlled current source: ``I(outp,outn) <+ gm*V(inp,inn)``."""
    instparams = [Parameter(name='gm', desc='Transconductance', unit='S',
                            default=1e-3)]

    @staticmethod
    def analog(inp, inn, outp, outn):
        bin_ = Branch(inp, inn)
        bout = Branch(outp, outn)
        return Contribution(bout.I, gm * bin_.V)


class VCVSHdl(Behavioural):
    """Voltage-controlled voltage source: ``V(outp,outn) <+ g*V(inp,inn)``."""
    instparams = [Parameter(name='g', desc='Voltage gain', unit='V/V',
                            default=1.0)]

    @staticmethod
    def analog(inp, inn, outp, outn):
        bin_ = Branch(inp, inn)
        bout = Branch(outp, outn)
        return Contribution(bout.V, g * bin_.V)


class DiodeHdl(Behavioural):
    """Ideal diode: ``I(p,n) <+ IS*(exp(V/vt) - 1)``.

    ``vt`` is an instance parameter (default the 300 K thermal voltage);
    ``elements.Diode`` reads temperature from ``epar`` instead -- epar
    access from HDL expressions is on the roadmap (hdl.md).  No limiting:
    Newton sees the raw exponential.
    """
    instparams = [Parameter(name='IS', desc='Saturation current', unit='A',
                            default=1e-13),
                  Parameter(name='vt', desc='Thermal voltage', unit='V',
                            default=0.02585199101160173)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I, IS * (sympy.exp(b.V / vt) - 1))


class VSinHdl(Behavioural):
    """Sinusoidal voltage source:
    ``V(p,n) <+ vo + va*sin(2*pi*freq*TIME + 2*pi*phase/360)``.

    Exercises the time-dependent source path: the ``TIME``-carrying terms
    compile into ``u(t)`` and their exact derivative into ``dudt(t)`` (the
    coupled-LTE path's requirement).
    """
    instparams = [Parameter(name='vo', desc='DC offset', unit='V',
                            default=0.0),
                  Parameter(name='va', desc='Amplitude', unit='V',
                            default=1.0),
                  Parameter(name='freq', desc='Frequency', unit='Hz',
                            default=1e3),
                  Parameter(name='phase', desc='Phase', unit='deg',
                            default=0.0)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        wt = 2 * sympy.pi * freq * TIME + 2 * sympy.pi * phase / 360
        return Contribution(b.V, vo + va * sympy.sin(wt))


class IdtHdl(Behavioural):
    """Integrator: ``V(op,on) <+ idt(V(ip,in))`` -- ``elements.Idt`` in
    one line.  The ``idt`` application compiles to a generated state
    equation; without an ``ic`` the DC solve is singular by LRM intent
    (use ``uic=True``), exactly like the hand-written element."""

    @staticmethod
    def analog(iplus, iminus, oplus, ominus):
        bi = Branch(iplus, iminus)
        bo = Branch(oplus, ominus)
        return Contribution(bo.V, idt(bi.V))


class IdtmodHdl(Behavioural):
    """Circular integrator: ``V(op,on) <+ idtmod(V(ip,in), ic, m, o)`` --
    ``elements.Idtmod`` in one line.

    The compiled element gets the whole idtmod treatment generated for
    free: the floored output wrap, the DC ic pin (fold-into-i, selected by
    ``epar.analysis_kind``), ``uic`` seeding of the wrapped ic, and the
    ``periodic_states()`` declaration that keeps the state bounded through
    the Phase-2 gauge shift (idtmod.md sections 5.1-5.2).
    """
    instparams = [Parameter(name='ic', desc='Initial value', unit='V',
                            default=0.0),
                  Parameter(name='modulus', desc='Output modulus',
                            unit='V/V', default=1.0),
                  Parameter(name='offset', desc='Offset', unit='V',
                            default=0.0)]

    @staticmethod
    def analog(iplus, iminus, oplus, ominus):
        bi = Branch(iplus, iminus)
        bo = Branch(oplus, ominus)
        return Contribution(bo.V, idtmod(bi.V, ic, modulus, offset))


## ======================================================================
## Self-heating -- roadmap item 1
## ======================================================================


class SelfHeating(object):
    """The self-heating thermal node, as a reusable block.

    An extra node whose "voltage" is the temperature rise ``dT`` above
    ambient in kelvin, an RC to the ambient reference, and the device's
    own dissipation injected into it as a current::

        I(th, tha) <+ V(th, tha)/rth + ddt(cth*V(th, tha))
        I(th, tha) <+ -p_dissipated

    KCL at ``th`` is then ``dT/rth + cth*d(dT)/dt = p``, so in steady
    state ``dT = p*rth`` and the thermal pole sits at ``1/(rth*cth)``.
    The device reads ``self.T`` -- ``$temperature + dT``, floored -- in
    place of ``TEMP`` and is otherwise unchanged.

    **A separate class per device, not a parameter on the existing one.**
    Three reasons, in order of force:

    1. *a parameter cannot switch the block off.*  An element compiles
       from **symbolic** parameters, so a `rth == 0` test is true at
       build time whatever the value; the thermal branch and its
       ``1/rth`` are in the expression regardless, and passing zero only
       evaluates them to infinity (`differentiable-numerics`: "a term
       switched OFF by a zero parameter is still evaluated").  What
       *does* switch it off is `Collapse`, which removes the branch
       before it is compiled -- and `Collapse` is declared in
       ``analog()``, i.e. it is part of the class, not of the instance;
    2. *the thermal node is a pin, and the pin count is the class.*
       Terminals come from the ``analog()`` signature and are fixed at
       class creation, so "with a thermal node" and "without" are
       necessarily two compiled elements;
    3. *making it a pin is worth more than making it internal.*  A
       terminal ``th`` can carry an external Foster or Cauer ladder built
       from ordinary ``R``/``C`` elements -- the roadmap's item 1b for
       free, no new operator -- and several devices can share one
       thermal network, which is the whole point of electrothermal
       simulation and is impossible with a node hidden inside one
       device.

    The ambient reference ``tha`` is a pin too, and has to be: an
    element can only reach the nodes in its own signature, so there is no
    way to say "referred to the global ground" from inside ``analog()``.
    Verilog-A gets this from the thermal discipline's own ground.
    Connect ``tha`` to ``gnd`` for an isothermal ambient.

    **One trap, and it is inherent to what `Collapse` can express.**
    `Collapse` may not merge two TERMINALS -- they belong to the parent
    circuit -- so ``rth = 0`` lands as a zero-volt source between ``th``
    and ``tha`` rather than as an absent branch.  For one device that is
    exactly right.  For several devices sharing one thermal node it is
    not: the ``rth = 0`` device's zero-volt source SHORTS the shared node
    to ambient and silently removes everybody else's temperature rise.
    A device that must not conduct heat needs its thermal pin left
    unconnected, not ``rth = 0``.
    ``test_two_devices_share_one_thermal_node`` pins the behaviour so
    nobody rediscovers it as a bug.

    Ten lines at the call site::

        heat = SelfHeating(th, tha, rth, cth)
        ...
        i = f(b.V, heat.T)
        return (Contribution(b.I, i),) + heat.dissipate(b.V*i)
    """

    def __init__(self, th, tha, rth, cth, tfloor=1.0):
        b = self.branch = Branch(th, tha)
        #: the temperature rise, in kelvin, as a named intermediate
        self.dT = _var(b.V, 'dT')
        ## FLOORED.  `dT` is a solution variable, so a Newton iterate can
        ## put the junction below absolute zero, and every temperature
        ## formula below takes `log(T/tnom)` and divides by `k*T/q`.  A
        ## smooth floor, not `Max`: `Max` differentiates to a step and
        ## the arms of everything downstream are evaluated anyway.
        ##
        ## The smoothing width is written RELATIVE to the floor, and it
        ## has to be small: `hypsmooth(x, e)` overshoots by `e^2/x`, so
        ## `e = 1e-3 K` biases a 300 K junction by 3.3e-9 K -- which
        ## sounds like nothing and is not, because the temperature enters
        ## through `exp(v/(n*k*T/q))` with an argument near 27, so it
        ## became a 2e-10 relative shift in the current and broke the
        ## claim that `rth = 0` reproduces the isothermal device exactly.
        ## At `1e-6*tfloor` the overshoot is below double roundoff at any
        ## ordinary temperature and the two agree bit for bit; the corner
        ## it rounds sits at 1 K, where only finiteness is being bought.
        self.T = _var(tfloor + _hypsmooth(TEMP + self.dT - tfloor,
                                          1e-6 * tfloor), 'Tj')
        self._network = (
            Contribution(b.I, b.V / rth + ddt(cth * b.V)),
            ## `rth = 0` is the isothermal limit, and it is the DEFAULT
            ## for a device that has not been characterised.  Collapsing
            ## the branch is what makes that expressible: the `1/rth`
            ## is never compiled for such an instance.  Both nodes are
            ## terminals here, so the collapse lands as a zero-volt
            ## source rather than a merged node (hdl.py's rule -- an
            ## element may not merge its parent's nodes), which pins
            ## `dT = 0` and costs one branch unknown.
            Collapse(b, rth <= 0))

    def dissipate(self, p):
        """The statements: the RC to ambient, the collapse rule, and the
        device's dissipated power ``p`` injected as a current into the
        thermal node."""
        return self._network + (Contribution(self.branch.I, -p),)


def _rtemp(r, tc1, tc2, dtnom):
    """``r*(1 + tc1*dT + tc2*dT^2)``, floored so it cannot reach zero.

    The polynomial is the standard SPICE/Verilog-A resistor tempco and it
    is not sign-definite: ``tc1 = -2e-3`` crosses zero at 500 K above
    nominal, and the model divides by it.  Floored at ``1e-6*r``, with a
    smoothing width written **relative to r** -- an absolute constant
    here would encode an assumption about a parameter we do not control
    (`differentiable-numerics`, "a smoothing constant is meaningless
    without its scale").
    """
    ## `floor + hypsmooth(x - floor, eps)` is the third hand-written
    ## clamp in this file -- the others are the power base in
    ## `_spice_diode` and the absolute-zero floor in `SelfHeating` -- and
    ## all three are `maxc(x, floor)` with a rounded corner.  `hdl.maxc`
    ## now exists; each of the three needed its own scan for `eps`, which
    ## is the part a helper cannot do for you.
    raw = _var(r * (1 + tc1 * dtnom + tc2 * dtnom ** 2), 'rraw')
    floor = 1e-6 * r
    return _var(floor + _hypsmooth(raw - floor, 1e-9 * r), 'rT')


def _thermal_params(rth=0.0, cth=0.0):
    return [Parameter(name='rth', desc='Thermal resistance junction to '
                      'ambient', unit='K/W', default=rth),
            Parameter(name='cth', desc='Thermal capacitance', unit='J/K',
                      default=cth)]


class RThermalHdl(Behavioural):
    """Self-heating resistor: ``I(p,n) <+ V(p,n)/r(T)`` with
    ``T = $temperature + V(th,tha)`` and ``V(th,tha)`` driven by the
    dissipation ``V(p,n)*I(p,n)``.

    Four terminals: the two electrical ones, then the thermal node and
    its ambient reference.  Tie ``th`` and ``tha`` together (or leave
    ``rth`` at its default of 0) for the isothermal resistor; hang an
    ``R``/``C`` ladder between them for a real package model.

    The fixed point it solves is ``dT = rth*V^2/r(1 + tc1*dT + ...)``,
    which for a positive ``tc1`` is self-limiting and for a negative one
    is the textbook thermal runaway -- both fall out of the Newton solve
    with an exact Jacobian, because the thermal node is a genuine
    unknown and not an outer iteration.
    """
    instparams = [Parameter(name='r', desc='Resistance at tnom',
                            unit='ohm', default=1e3),
                  Parameter(name='tc1', desc='Linear tempco', unit='1/K',
                            default=0.0),
                  Parameter(name='tc2', desc='Quadratic tempco',
                            unit='1/K^2', default=0.0),
                  Parameter(name='tnom', desc='Parameter measurement '
                            'temperature', unit='K', default=300.15)] \
        + _thermal_params()

    @staticmethod
    def analog(plus, minus, th, tha):
        heat = SelfHeating(th, tha, rth, cth)                      # noqa
        b = Branch(plus, minus)
        rT = _rtemp(r, tc1, tc2, _var(heat.T - tnom, 'dtnom'))     # noqa
        i = _var(b.V / rT, 'i')
        return (Contribution(b.I, i),) + heat.dissipate(b.V * i)


## ======================================================================
## The SPICE level-1 junction diode -- roadmap item 2
## ======================================================================


def _spice_diode_params():
    """The SPICE diode card, in SPICE's own names and defaults.

    ``is`` is a Python keyword, so the saturation current keeps
    ``DiodeHdl``'s spelling ``IS``.
    """
    return [
        Parameter(name='IS', desc='Saturation current at tnom', unit='A',
                  default=1e-14),
        Parameter(name='rs', desc='Ohmic series resistance', unit='ohm',
                  default=0.0),
        Parameter(name='n', desc='Emission coefficient', unit='',
                  default=1.0),
        Parameter(name='tt', desc='Transit time', unit='s', default=0.0),
        Parameter(name='cjo', desc='Zero-bias junction capacitance',
                  unit='F', default=0.0),
        Parameter(name='vj', desc='Junction potential at tnom', unit='V',
                  default=1.0),
        Parameter(name='m', desc='Grading coefficient', unit='',
                  default=0.5),
        Parameter(name='eg', desc='Activation energy', unit='eV',
                  default=1.11),
        Parameter(name='xti', desc='Saturation-current temperature '
                  'exponent', unit='', default=3.0),
        Parameter(name='fc', desc='Forward-bias depletion coefficient',
                  unit='', default=0.5),
        ## SPICE's BV default is "no breakdown".  1e30 V is that, and it
        ## is the arithmetic spelling of it: the breakdown exponent's
        ## argument goes to -1e34 and its `exp` underflows to exactly
        ## 0.0, which is the one benign IEEE condition (hdl.md 3.2c).
        Parameter(name='bv', desc='Reverse breakdown voltage', unit='V',
                  default=1e30),
        Parameter(name='ibv', desc='Current at breakdown', unit='A',
                  default=1e-3),
        Parameter(name='kf', desc='Flicker-noise coefficient', unit='',
                  default=0.0),
        Parameter(name='af', desc='Flicker-noise exponent', unit='',
                  default=1.0),
        Parameter(name='area', desc='Area scaling factor', unit='',
                  default=1.0),
        Parameter(name='tnom', desc='Parameter measurement temperature',
                  unit='K', default=300.15),
    ]


def _spice_diode(a, c, T, IS, rs, n, tt, cjo, vj, m, eg, xti, fc, bv, ibv,
                 kf, af, area, tnom):
    """The junction itself, shared by the isothermal and self-heating
    classes.  Returns ``(statements, p_dissipated)``.

    Every parameter is passed **explicitly** rather than read from the
    module globals the compiler injects them into, because those globals
    belong to whichever class was compiled last (see the import block).

    Equations: Massobrio & Antognetti, *Semiconductor Device Modeling
    with SPICE*, ch. 1 -- the same set ngspice's `diode` implements.

    * ``Eg(T) = 1.16 - 7.02e-4*T^2/(T + 1108)``
    * ``IS(T) = IS*exp((T/Tnom - 1)*eg/(n*Vt) + (xti/n)*ln(T/Tnom))``
    * ``VJ(T) = VJ*T/Tnom - 3*Vt*ln(T/Tnom) - Eg(Tnom)*T/Tnom + Eg(T)``
    * ``CJ0(T) = CJ0*(1 + m*(4e-4*(T - Tnom) - (VJ(T)/VJ - 1)))``
    * depletion charge with SPICE's ``fc`` linearisation above
      ``fc*VJ(T)``, C1 at the seam;
    * diffusion charge ``tt*I_fwd``;
    * breakdown as an additive ``-ibv*exp(-(v + bv)/Vt)``, which is
      ``-ibv`` exactly at ``v = -bv`` and smooth everywhere -- SPICE's
      own three-region form is not;
    * shot noise ``2*q*|I|`` and flicker noise ``kf*|I|^af/f`` on the
      junction, thermal noise ``4*k*T/rs`` on the series resistance.
    """
    di = Node('di')
    brs = Branch(a, di)
    bd = Branch(di, c)

    ## -- temperature path ----------------------------------------------
    vtT = _var(_vt(T), 'vtT')
    tr = _var(T / tnom, 'tratio')
    ltr = _var(sympy.log(tr), 'ltratio')
    egT = _var(1.16 - 7.02e-4 * T ** 2 / (T + 1108.0), 'egT')
    egn = _var(1.16 - 7.02e-4 * tnom ** 2 / (tnom + 1108.0), 'egtnom')
    ## `expl`, not `exp`.  The argument is bounded ABOVE by
    ## `eg*q/(n*k*tnom)` (~43 for a silicon card at n=1), so on any
    ## sensible card the two are identical; `expl` costs nothing and
    ## keeps a card with n < 0.5 -- or a runaway thermal node -- finite.
    isT = _var(IS * _expl((tr - 1) * eg / (n * vtT) + xti / n * ltr),
               'isT')
    vjT = _var(vj * tr - 3 * vtT * ltr - egn * tr + egT, 'vjT')
    cjT = _var(cjo * (1 + m * (4e-4 * (T - tnom) - (vjT / vj - 1))), 'cjT')

    ## -- static current ------------------------------------------------
    vd = _var(bd.V, 'vd')
    nvt = _var(n * vtT, 'nvt')
    ifwd = _var(isT * (_expl(vd / nvt) - 1), 'ifwd')
    ibrk = _var(ibv * _expl(-(vd + bv) / vtT), 'ibrk')
    idio = _var(area * (ifwd - ibrk), 'id')

    ## -- charge --------------------------------------------------------
    fcvj = _var(fc * vjT, 'fcvj')
    ubar = _var(1 - vd / vjT, 'ubar')
    ## FLOOR THE BASE OF THE POWER.  `ubar**(1-m)` is NaN for `vd > vjT`,
    ## and the low arm is evaluated at every bias whether it is selected
    ## or not.  The floor is placed at HALF the seam value `1-fc`, not at
    ## the seam and not at zero: put it at the seam and the smoothing
    ## perturbs the arm exactly where C1 continuity is being claimed;
    ## put it at zero and the derivative of a fractional power diverges
    ## there.  `ubar = ufl` is `vd = vjT*(1+fc)/2`, which for any
    ## `fc < 1` is ABOVE the seam `fc*vjT` -- so the corner provably sits
    ## in the arm that is always discarded, and its sharpness costs
    ## nothing.
    ##
    ## The smoothing width is relative to `1-fc`, and it perturbs `ubase`
    ## at the seam by about `c^2`: scanned, `c = 1e-6` left a 6e-25 C
    ## offset on a charge that is exactly zero at zero bias -- physically
    ## irrelevant, but visible against an independent reference and
    ## therefore a real bias.  `c = 1e-9` puts it at 1e-18, below double
    ## roundoff, and the corner it rounds is unreachable anyway.
    ##
    ## THIS IS THE WORKAROUND THAT `hdl.safe_pow` REPLACES.  Three lines
    ## and a floor-placement argument, written by hand because the
    ## kernel had no floored power when this model was written; the
    ## whole block is `safe_pow(1 - vd/vjT, 1 - m, lo=0.5*(1 - fc))`, and
    ## `safe_pow`'s hard `maxc` clamp is the better instrument here
    ## because the floor sits strictly inside the domain.  Kept as
    ## written so that the friction it records stays visible.
    ufl = _var(0.5 * (1 - fc), 'ufloor')
    ubase = _var(ufl + _hypsmooth(ubar - ufl, 1e-9 * (1 - fc)), 'ubase')
    qlo = _var(cjT * vjT * (1 - ubase ** (1 - m)) / (1 - m), 'qdeplo')
    f1 = _var(vjT * (1 - (1 - fc) ** (1 - m)) / (1 - m), 'f1')
    f2 = _var((1 - fc) ** (1 + m), 'f2')
    f3 = _var(1 - fc * (1 + m), 'f3')
    qhi = _var(cjT * (f1 + (f3 * (vd - fcvj)
                            + m / (2 * vjT) * (vd ** 2 - fcvj ** 2)) / f2),
               'qdephi')
    qdep = _var(sympy.Piecewise((qlo, vd < fcvj), (qhi, True)), 'qdep')
    ## Diffusion charge from the FORWARD current only: the breakdown
    ## term is not minority-carrier injection and carries no transit
    ## time, and `tt*ibrk` at 1 kV reverse would be a large spurious
    ## capacitance.
    qtot = _var(area * (qdep + tt * ifwd), 'qtot')

    ## -- noise ---------------------------------------------------------
    ## `safe_abs`, not `Abs`: the PSD must be non-negative and the
    ## current is negative over most of the reverse range.  (`CY` is
    ## evaluated, never differentiated, so the regularisation here is
    ## about the value only.)  With `kf = 0` the flicker term is zero and
    ## still evaluated -- which is fine, because `|I|^af` is finite for
    ## every bias, floored below by `safe_abs`'s own eps.
    iabs = _var(_safe_abs(idio), 'iabs')
    noise = (Contribution(bd.I, _white_noise(2 * _QE * iabs)),
             Contribution(bd.I, _flicker_noise(kf * iabs ** af, 1)),
             ## Thermal noise of the series resistance.  It goes away
             ## with the branch when rs collapses, which is the correct
             ## behaviour and comes for free.
             Contribution(brs.I, _white_noise(4 * _KB * T * area / rs)))

    ## -- statements ----------------------------------------------------
    stmts = noise + (
        Contribution(brs.I, brs.V * area / rs),
        ## The Verilog-A optional-parasitic idiom, and the reason
        ## `rs = 0` (SPICE's default) is not a division by zero: the
        ## collapsed variant never compiles the `1/rs`.
        Collapse(brs, rs <= 0),
        Contribution(bd.I, idio),
        Contribution(bd.I, ddt(qtot)))
    ## Total static dissipation: the series chain carries one current, so
    ## the terminal-to-terminal voltage times it is the whole of it.  The
    ## charge term is deliberately absent -- storage does not dissipate.
    return stmts, _var(Branch(a, c).V * idio, 'pdiss')


class DiodeSpiceHdl(Behavioural):
    """The full SPICE level-1 junction diode.

    ``DiodeHdl`` above is the bare exponential; this is the card every
    PDK actually ships -- series resistance on an internal node,
    depletion and diffusion charge, the ``is``/``vj``/``cj0``
    temperature path, reverse breakdown and area scaling -- and it is
    the model that exercises ``ddt``, `Collapse`, the internal-node
    machinery and the temperature path in one element.

    Parameters are SPICE's, with SPICE's defaults (``rs = 0``,
    ``cjo = 0``, ``tt = 0``, no breakdown), so a bare
    ``DiodeSpiceHdl('a', 'b')`` is the ideal exponential and each
    optional block is switched on by giving its parameter.  ``rs = 0``
    collapses its branch away rather than dividing by zero;
    ``cjo = tt = 0`` really does leave the charge at zero, because the
    charge expression is a product and not a quotient.

    Shot and flicker noise on the junction and thermal noise on ``rs``
    are included (roadmap item 14, for this device).  There is no
    ``$limit``, and no PCNR participation either -- the layer refuses a
    device that carries charge -- so the ``expl`` in the forward term is
    the whole of what keeps a wild Newton iterate finite.  See the
    import block for why it is ``expl`` and not ``limexp``.
    """
    instparams = _spice_diode_params()

    @staticmethod
    def analog(plus, minus):
        return _spice_diode(plus, minus, TEMP, IS, rs, n, tt, cjo, vj,  # noqa
                            m, eg, xti, fc, bv, ibv, kf, af,            # noqa
                            area, tnom)[0]                              # noqa


class DiodeSpiceThermalHdl(Behavioural):
    """`DiodeSpiceHdl` with the thermal node attached.

    The diff against `DiodeSpiceHdl.analog` is four lines: build the
    `SelfHeating` block, pass its ``T`` instead of ``TEMP``, and append
    its statements.  That is the roadmap's claim for item 1 -- "~10
    lines per device, and it retroactively upgrades every model in the
    library" -- measured on the second model in the library.

    Terminals ``(plus, minus, th, tha)``.  With ``rth`` at its default
    of 0 the thermal branch collapses to a zero-volt source and this
    element is `DiodeSpiceHdl` again, to the last digit.
    """
    instparams = _spice_diode_params() + _thermal_params()

    @staticmethod
    def analog(plus, minus, th, tha):
        heat = SelfHeating(th, tha, rth, cth)                          # noqa
        stmts, p = _spice_diode(plus, minus, heat.T, IS, rs, n, tt,    # noqa
                                cjo, vj, m, eg, xti, fc, bv, ibv,      # noqa
                                kf, af, area, tnom)                    # noqa
        return stmts + heat.dissipate(p)
