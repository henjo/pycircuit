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
``XtalHdl``      Butterworth-Van Dyke quartz resonator: a series R-L
                 written as a potential contribution, two ``ddt``
                 branches, and the two closed-form resonances an
                 oscillator is built around
``FerriteBeadHdl``  ``laplace_nd`` from network coefficients -- the
                 parallel R-L-C every bead datasheet is a fit to
``RSkinHdl``     ``laplace_zp``: an interlaced pole/zero ladder whose
                 average slope is HALF, i.e. a rational approximation
                 to the skin effect's ``sqrt(s)``
``ComparatorHdl``  the first production user of ``Cross``: hysteresis
                 as positive feedback (the latch is in the *solution*,
                 not in a stored state) with the crossing declared at
                 the closed-form fold voltage
``MemristorHdl``  ``idt`` as a real DAE state -- and the model that
                 shows what is missing, since a state that appears in
                 its own derivative cannot be written with ``idt`` at
                 all
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
                                   Collapse, Cross, ddt, idt, idtmod,
                                   laplace_nd, laplace_zp, TIME, TEMP)

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
from pycircuit.circuit.hdl import safe_div as _safe_div
from pycircuit.circuit.hdl import maxc as _maxc
from pycircuit.circuit.hdl import minc as _minc
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


## ======================================================================
## Passive resonators and dispersive impedances -- roadmap items 5 and 8
## ======================================================================


def bvd_from_fs_q(fs, q, cm, c0):
    """The four BVD elements from what a crystal datasheet actually prints.

    A datasheet gives the series resonance ``fs``, the quality factor
    ``q`` (or the ESR), the motional capacitance ``cm`` (often called
    ``C1``, of order femtofarads) and the shunt capacitance ``c0``.  It
    never gives the motional inductance, because it is millihenries and
    nobody would believe it.  The two identities are

    * ``lm = 1/((2*pi*fs)^2 * cm)`` -- the series resonance;
    * ``rm = 2*pi*fs*lm/q = 1/(2*pi*fs*cm*q)`` -- the definition of Q for
      a series RLC.

    Returns the keyword dict `XtalHdl` takes.
    """
    ws = 2 * 3.141592653589793 * fs
    lm = 1.0 / (ws * ws * cm)
    return dict(lm=lm, cm=cm, rm=ws * lm / q, c0=c0)


class XtalHdl(Behavioural):
    """Quartz crystal / ceramic resonator -- the Butterworth-Van Dyke model.

    The one element in the repo that can set an oscillator's frequency.
    A motional arm ``rm``-``lm``-``cm`` in series, shunted by the holder
    capacitance ``c0``::

        p o--+--[ rm ]--[ lm ]--(m)--[ cm ]--+--o n
             |                               |
             +------------[ c0 ]-------------+

    Written as three contributions: the series R-L is one branch whose
    **potential** is contributed (``V = rm*I + ddt(lm*I)``, the inductor
    idiom from `LHdl` with a resistive term added), and ``cm`` and ``c0``
    are ordinary ``ddt`` current contributions.  The internal node ``m``
    is the junction between the inductor and the motional capacitor.

    The two resonances are the reason the model exists, and both are
    closed form:

    * **series** ``fs = 1/(2*pi*sqrt(lm*cm))`` -- the motional arm's
      reactance vanishes and the impedance falls to ``rm``;
    * **parallel** ``fp = fs*sqrt(1 + cm/c0)`` -- the motional arm is
      inductive and cancels ``c0``; the impedance peaks at roughly
      ``lm/(rm*(c0 + cm))``.

    The *pulling range* ``fp - fs`` is what an oscillator's load
    capacitance moves the frequency over, and it is set entirely by the
    capacitance ratio ``c0/cm`` -- 200 to 400 for quartz, which is why a
    crystal pulls a few hundred ppm and an LC tank does not.

    Q is ``2*pi*fs*lm/rm``; for the default card it is 1.3e5, which is
    also a warning: a transient that has to resolve a 1e5-Q resonance
    needs to run for 1e5 cycles before the envelope settles.  AC is the
    analysis this element is for.  `bvd_from_fs_q` converts a datasheet.

    No parameter is ever a divisor here, so every card is finite --
    including ``rm = 0`` (a lossless resonator) and ``c0 = 0`` (a bare
    series arm).  That is not luck: it is what writing the motional arm
    as a *potential* contribution buys, and it is why there is no
    `Collapse` in this model.
    """
    instparams = [Parameter(name='lm', desc='Motional inductance', unit='H',
                            default=8.0e-3),
                  Parameter(name='cm', desc='Motional capacitance', unit='F',
                            default=3.2e-15),
                  Parameter(name='rm', desc='Motional resistance (ESR)',
                            unit='ohm', default=12.0),
                  Parameter(name='c0', desc='Shunt (holder) capacitance',
                            unit='F', default=3.5e-12)]

    @staticmethod
    def analog(plus, minus):
        mid = Node('m')
        bmot = Branch(plus, mid, 'mot')
        bcm = Branch(mid, minus, 'cm')
        bc0 = Branch(plus, minus, 'c0')
        return (Contribution(bmot.V, rm * bmot.I + ddt(lm * bmot.I)),  # noqa
                Contribution(bcm.I, ddt(cm * bcm.V)),                  # noqa
                Contribution(bc0.I, ddt(c0 * bc0.V)))                  # noqa


class FerriteBeadHdl(Behavioural):
    """Ferrite bead / lossy inductor, as the parallel R-L-C a datasheet fits.

    A bead is an inductor below self-resonance, a resistor at it, and a
    capacitor above it -- the impedance curve every bead datasheet plots.
    That is exactly a parallel ``rp || ls || cp`` in series with the dc
    winding resistance::

        Z(s) = rdc + s*ls / (1 + s*ls/rp + s^2*ls*cp)

    which is written here as one `laplace_nd` on the branch **current**::

        Contribution(b.V, rdc*b.I
                     + laplace_nd(b.I, [0, ls], [1, ls/rp, ls*cp]))

    `laplace_nd` rather than `laplace_zp` and the reason is worth
    recording: the poles are

        s = (-1 +- sqrt(1 - 4*rp^2*cp/ls)) / (2*rp*cp)

    a square root of a parameter expression whose RADICAND changes sign
    -- real poles for a low-Q bead (the default card, ``Q = 0.32``), a
    conjugate pair for a high-Q one.  Pole-zero form would need that
    radical symbolic in every instance, with a branch on its sign;
    coefficient form is three products.  **Whenever
    the poles come from a physical network rather than from a fit,
    `laplace_nd` is the honest spelling** -- `laplace_zp` earns its keep
    when the roots are what you know (see `RSkinHdl`).

    The three landmarks, all closed form and all independent of the
    realisation:

    * ``|Z| -> 2*pi*f*ls`` well below self-resonance;
    * ``f0 = 1/(2*pi*sqrt(ls*cp))``, where ``|Z| = rdc + rp`` exactly and
      the phase is zero -- this is the "impedance at 100 MHz" a
      datasheet quotes;
    * ``|Z| -> 1/(2*pi*f*cp)`` above it, the reason a bead stops working
      where you most want it to.

    Defaults are a 600 ohm-at-250-MHz signal bead.  ``rdc`` is a real
    term, not a numerical guard: it is what limits the dc voltage drop
    at the supply current, and it is the only thing the model has at
    ``s = 0`` (the parallel L shorts the rest).
    """
    instparams = [Parameter(name='ls', desc='Bead inductance', unit='H',
                            default=1.2e-6),
                  Parameter(name='rp', desc='Impedance at self-resonance',
                            unit='ohm', default=600.0),
                  Parameter(name='cp', desc='Parallel (winding) capacitance',
                            unit='F', default=0.35e-12),
                  Parameter(name='rdc', desc='DC winding resistance',
                            unit='ohm', default=0.08)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        ## Ascending powers of s, as the LRM specifies for laplace_nd.
        ## The numerator's leading 0 is the zero at the origin: a bead
        ## is a SHORT at dc apart from `rdc`, and writing it as a
        ## coefficient rather than as a `laplace_zp` root at zero keeps
        ## the transfer function dimensionless-times-henries instead of
        ## silently acquiring a 1/time (see `laplace_zp`'s convention).
        return Contribution(b.V, rdc * b.I                             # noqa
                            + laplace_nd(b.I, [0, ls],                 # noqa
                                         [1, ls / rp, ls * cp]))       # noqa


#: Pole/zero pairs in the default `RSkinHdl` ladder, and the decades of
#: normalised frequency the fit is placed over.  STRUCTURAL, not
#: parameters: `laplace_zp` takes a python LIST of roots, so its length
#: is fixed when `analog()` runs, i.e. once per class -- the same reason
#: `elements.VPWL` cannot move to the DSL (roadmap sec. 5).  Use
#: `skin_effect_resistor()` for a different order.
_SKIN_PAIRS, _SKIN_LO, _SKIN_HI = 11, -5.0, 5.0


def _fractional_zp(alpha, lo, hi, pairs):
    """Oustaloup's recursive approximation to ``s**alpha``, normalised.

    Returns ``(zeros, poles, gain)`` for a rational ``H`` with
    ``gain*H(s) ~ s**alpha`` over ``10**lo < |s| < 10**hi``, with
    ``zeros``/``poles`` in `laplace_zp`'s flat ``[re, im, re, im, ...]``
    form and ``H(0) = 1`` (`laplace_zp`'s normalisation, which is not
    Oustaloup's -- hence the returned ``gain``).

    Oustaloup, Levron, Mathieu & Nanot, "Frequency-band complex
    noninteger differentiator: characterization and synthesis", *IEEE
    Trans. Circuits Syst. I* **47**(1):25-39, 2000, eq. (18):

        w_z[k] = wb*(wh/wb)**((k + n + (1-alpha)/2)/(2n+1))
        w_p[k] = wb*(wh/wb)**((k + n + (1+alpha)/2)/(2n+1))
        K      = wh**alpha,        k = -n .. n

    i.e. ``2n+1`` real poles and zeros interlaced geometrically, the
    zeros displaced below the poles by ``alpha`` of a cell.  A cell's
    zero-pole pair contributes ``alpha`` of a decade of gain per decade,
    which is the whole trick: **a fractional slope is an interlaced
    ladder, and the fraction is the interlacing offset.**

    ACCURACY, measured (`test_elements_hdl_library2`): with the defaults
    -- 11 pairs over 10 decades -- the magnitude is within 1.7% and the
    phase within 1.4 degrees of ``s**0.5`` over the middle 8 decades.
    The outermost half-decade at each end is where the ladder runs out
    and must not be used; that is what ``lo``/``hi`` being *wider* than
    the advertised band is for.
    """
    import numpy as _np
    n = (pairs - 1) // 2
    ks = _np.arange(-n, n + 1)
    wb, wh = 10.0 ** lo, 10.0 ** hi
    ratio = wh / wb
    wz = wb * ratio ** ((ks + n + (1 - alpha) / 2) / pairs)
    wp = wb * ratio ** ((ks + n + (1 + alpha) / 2) / pairs)
    ## Oustaloup's K is for the form `prod (s + wz)/(s + wp)`;
    ## `laplace_zp` builds `prod (1 - s/root)`, which is the same
    ## function divided by `prod wz/wp`.  Convert once, here, rather
    ## than leaving a factor for every caller to rediscover.
    gain = float(wh ** alpha * _np.prod(wz / wp))
    zeros, poles = [], []
    for z, p in zip(wz, wp):
        ## A REAL root, so the imaginary member of the pair is 0.0 and
        ## the root itself is NEGATIVE: `laplace_zp`'s factor is
        ## `1 - s/root`, so a stable pole is a root in the left half
        ## plane and writing `+wp` here would silently build an
        ## unstable filter that still compiles and still solves at dc.
        zeros += [-float(z), 0.0]
        poles += [-float(p), 0.0]
    return zeros, poles, gain


def skin_effect_resistor(pairs=_SKIN_PAIRS, lo=_SKIN_LO, hi=_SKIN_HI):
    """Build a skin-effect resistor class with a ladder of ``pairs``
    pole/zero pairs placed over ``10**lo .. 10**hi`` times ``2*pi*fs``.

    The order is a property of the CLASS, not of an instance -- see
    `_SKIN_PAIRS`.  `RSkinHdl` is ``skin_effect_resistor()``.
    """
    zeros, poles, gain = _fractional_zp(0.5, lo, hi, pairs)

    class RSkin(Behavioural):
        __doc__ = """Skin-effect resistor: ``Z(s) = rdc*(1 + sqrt(s/ws))``.

    Above the onset frequency the current crowds into a skin one
    penetration depth thick, so the resistance rises as ``sqrt(f)`` and
    an equal internal reactance rises with it -- the impedance phase
    tends to **45 degrees**, not to 90.  That half-power law is not a
    rational function of ``s``, and this element is the first thing in
    the repo that needs `laplace_zp` for what it is for: an interlaced
    ladder of %d real pole/zero pairs whose average slope is half.

    ``ws = 2*pi*fs``, and ``fs`` is defined as **the frequency at which
    the skin term equals the dc resistance** -- ``|sqrt(j*w/ws)| = 1``
    there, so ``|Z(fs)| = rdc*|1 + sqrt(j)| = 1.8478*rdc`` and the phase
    is 22.5 degrees.  For a round wire of radius ``a`` in a material of
    conductivity ``sigma``, ``fs`` is where the penetration depth
    ``sqrt(2/(w*mu*sigma))`` is about ``a/2``.

    BAND.  The ladder is fitted over ``fs*1e%+.0f`` to ``fs*1e%+.0f``
    and is accurate over the middle of that: measured, better than 2%%
    in magnitude and 1.5 degrees in phase from ``fs*1e-4`` to
    ``fs*1e+3``.  Outside it the element is still a well-behaved
    passive impedance -- it flattens to a resistor at each end rather
    than doing anything exciting -- but it is no longer this physics.
    **A rational approximation of a fractional power has a band, and
    the band is part of the model**; ask for a different one with
    `skin_effect_resistor`.

    The one artefact worth knowing: ``H(0) = 1`` for a pole-zero
    product, and ``sqrt(0) = 0``, so the ladder cannot make the skin
    term vanish at dc.  It leaves ``gain*H(0) = %.3e`` of it, i.e. the
    dc resistance is ``rdc*(1 + %.3e)`` and not ``rdc``.  That is
    0.%03d%% on the default ladder, it is exactly computable rather
    than approximately small, and `test_skin_dc_floor_is_the_ladders_own`
    pins it.
    """ % (pairs, lo, hi, gain, gain, int(round(gain * 1e5)))

        instparams = [Parameter(name='rdc', desc='DC resistance', unit='ohm',
                                default=1.0),
                      Parameter(name='fs', desc='Frequency at which the skin '
                                'term equals rdc', unit='Hz', default=1e6)]

        @staticmethod
        def analog(plus, minus):
            b = Branch(plus, minus)
            ## The ladder is computed once, in NORMALISED frequency, and
            ## scaled by `ws` here.  Two reasons, and the second is the
            ## one that decides it:
            ##
            ## 1. `H(s/ws)` is the same filter for every instance, so
            ##    the pole positions are `number * ws` and sympy's
            ##    expansion of the order-11 product collects into
            ##    coefficients `c_k/ws**k` with NUMERIC `c_k`.  Compile
            ##    time is 0.3 s;
            ## 2. with the band exponents symbolic too, the same product
            ##    expands into 2**11 symbolic terms and the class does
            ##    not finish compiling.  The SHAPE of a laplace ladder
            ##    has to be numeric; only its scale can be a parameter.
            ws = 2 * sympy.pi * fs                                     # noqa
            zs = [z * ws if k % 2 == 0 else z
                  for k, z in enumerate(zeros)]
            ps = [p * ws if k % 2 == 0 else p
                  for k, p in enumerate(poles)]
            return Contribution(b.V, rdc * (b.I + gain                 # noqa
                                            * laplace_zp(b.I, zs, ps)))

    RSkin.__name__ = RSkin.__qualname__ = 'RSkin%d' % pairs
    return RSkin


RSkinHdl = skin_effect_resistor()
RSkinHdl.__name__ = RSkinHdl.__qualname__ = 'RSkinHdl'


## ======================================================================
## Behavioural: the comparator with hysteresis -- roadmap item 7
## ======================================================================


class ComparatorHdl(Behavioural):
    """Comparator with hysteresis, and the repo's first user of `Cross`.

    Four terminals ``(inp, inn, outp, outn)``::

        V(outp,outn) <+ vol + (voh - vol)/2 * (1 + q)
        q            =  tanh(gain*(V(inp,inn) - vref + vhyst*q))

    **The hysteresis is positive feedback, not stored state**, and that
    is the whole design.  ``q`` -- the normalised decision, in (-1,1) --
    appears on both sides of its own equation, so for a loop gain
    ``gain*vhyst > 1`` the equation has three roots over a window of
    input voltage and Newton stays on the branch it is already on.  A
    real Schmitt trigger latches for exactly this reason, and it means
    the element needs no event, no discrete state and no memory the
    solver cannot see: the latch is in the *solution*, and rejecting a
    timestep un-latches it correctly for free.

    THE THRESHOLDS ARE NOT ``vref +- vhyst``.  They are the fold points
    of that cubic-like curve, where ``dq/dvin`` becomes infinite:

        qf   = sqrt(1 - 1/(gain*vhyst))
        vth  = vhyst*qf - atanh(qf)/gain

    and the loop is ``vref +- vth``, symmetric.  ``vth -> vhyst`` only
    as ``gain -> infinity``; at ``gain*vhyst = 4`` it is 0.537 of it, so
    a model that advertised ``vhyst`` as the trip point would be wrong
    by a factor of two.  For ``gain*vhyst <= 1`` there is no fold, the
    element is a plain saturating comparator, and the formula degrades
    to ``vth = 0`` on its own -- which is the right answer.

    ``Cross`` is declared at those two thresholds, one per direction,
    and this is the first place in the repo where the difference is
    measurable.  Measured (`test_cross_lands_the_comparator_edges`) on a
    1 kHz sine into the default card: the nearest timepoint to a
    threshold crossing goes from 7.4e-7..1.2e-5 s away to **3e-16 s on
    three edges and 6.6e-7 on the fourth**.

    That scatter is the finding, not a blemish.  `Cross` extrapolates
    linearly from the last two accepted points and publishes ONE
    prediction; a second, far better one only arrives if the controller
    polls again after landing near the corner, and there is no
    iteration to convergence.  So the worst edge of a run is the one
    the predictor got only one shot at, and 6.6e-7 s is exactly what a
    first-order prediction over a 4e-5 s step on a sine is worth.  A
    caller should read `Cross` as "brackets the corner to first order",
    which is what its docstring says, and not as "lands on it" -- even
    though three times out of four it did.

    It costs about **3x the accepted steps**, because every breakpoint
    restarts the step controller at a small step and it has to grow
    back; the existing `test_cross_lands_timepoints_on_the_crossing`
    sees no such cost because its comparator does not latch and so has
    no corner to resolve.  All three numbers are the honest report of
    what `Cross` is for.

    The model is written with `var()`, i.e. on the **chained** compile
    path, and that is a deliberate regression guard rather than an
    optimisation -- see the comment in ``analog()``.

    ``q`` lives on an internal node ``lat``, referred to ``outn``, put
    there by a **potential** contribution.  Nothing else touches that
    node, so KCL there forces its branch current to zero and the node
    injects exactly nothing into ``outn``: it is a scratch signal, not a
    circuit.  That idiom -- an internal node written by a V-contribution
    and read back -- is how a Verilog-A local variable that must be
    *shared with the solver* is spelled here, and `MemristorHdl` needs
    it for a different reason.
    """
    instparams = [Parameter(name='vref', desc='Reference (trip centre)',
                            unit='V', default=0.0),
                  Parameter(name='vhyst', desc='Half hysteresis at infinite '
                            'gain', unit='V', default=0.2),
                  Parameter(name='gain', desc='Small-signal gain', unit='1/V',
                            default=20.0),
                  Parameter(name='voh', desc='Output high', unit='V',
                            default=1.0),
                  Parameter(name='vol', desc='Output low', unit='V',
                            default=-1.0)]

    @staticmethod
    def analog(inp, inn, outp, outn):
        bin_ = Branch(inp, inn)
        bout = Branch(outp, outn)
        lat = Node('lat')
        blat = Branch(lat, outn, 'q')
        q = blat.V
        ## The fold voltage, in closed form.  Both guards earn their
        ## place and neither is a smoothing:
        ##
        ## * `safe_div(., lg, 1e-10)` because `vhyst = 0` (no hysteresis)
        ##   is a legitimate card and `1/lg` would be a divide-by-zero
        ##   raised on EVERY evaluation of the crossing function, not
        ##   just on that card -- both arms of everything are evaluated.
        ##   The eps is relative to a DIMENSIONLESS loop gain, which is
        ##   the one place an absolute constant is defensible;
        ## * `maxc(lg - 1, 0)` because below unity loop gain the radicand
        ##   is negative and `sqrt` of it is NaN.  A hard clamp, not a
        ##   smoothing: the floor sits at the fold's own disappearance,
        ##   the physics below it really is "no fold", and the derivative
        ##   of the clamped arm is an honest zero rather than the `0/0`
        ##   a `sqrt(max(x,0))` would give (`differentiable-numerics`).
        lg = _var(gain * vhyst, 'loopgain')                            # noqa
        qf = _var(sympy.sqrt(_safe_div(_maxc(lg - 1, 0.0), lg, 1e-10)),
                  'qfold')
        vth = _var(vhyst * qf - sympy.atanh(qf) / gain, 'vth')         # noqa
        ## `var()`, i.e. this model is deliberately CHAINED, and that is
        ## not an optimisation -- three intermediates would inline
        ## perfectly well.  It is a regression guard.  `Cross` was
        ## installed on the class AFTER an early `return` taken by every
        ## model that calls `var()`, so from the day `@cross` was built
        ## until 2026-08-25 it silently did nothing on any model a real
        ## author would write -- while `explain()` went on reporting
        ## "1 @cross event", because it reads the compiler's record
        ## rather than the class.  A capability with no production user
        ## is a capability nobody has run.  This one now has one, on the
        ## path where it was broken.
        return (Contribution(blat.V,
                             sympy.tanh(_var(gain * (bin_.V - vref     # noqa
                                                     + vhyst * q),     # noqa
                                             'arg'))),
                Contribution(bout.V,
                             vol + (voh - vol) / 2 * (1 + q)),         # noqa
                ## One crossing per threshold, each with its own
                ## direction: the rising edge can only happen at the
                ## upper fold and the falling edge only at the lower, so
                ## a direction-0 pair would ask the controller to land on
                ## two corners that cannot occur.
                Cross(bin_.V - vref - vth, +1),                        # noqa
                Cross(bin_.V - vref + vth, -1))


## ======================================================================
## The memristor -- roadmap item 11
## ======================================================================


class MemristorHdl(Behavioural):
    """HP / Strukov memristor with the Biolek window -- ``idt`` as a
    genuine DAE state.

    Strukov, Snider, Stewart & Williams, "The missing memristor found",
    *Nature* **453**:80-83, 2008, with the boundary window of Biolek,
    Biolek & Biolkova, "SPICE model of memristor with nonlinear dopant
    drift", *Radioengineering* **18**(2):210-214, 2009::

        R(x) = ron*x + roff*(1 - x)                 x in [0, 1]
        i    = V(p,n)/R(x)
        dx/dt = muv*ron/d^2 * i * f(x, i)
        f(x, i) = 1 - (x - stp(-i))^(2p)            (Biolek)

    ``x`` is the doped fraction of the film.  The signature is a
    **pinched hysteresis loop**: the i-v curve passes exactly through
    the origin (no current without voltage -- it is a resistor, not a
    source) and opens into two lobes whose area shrinks as frequency
    rises, because the state has less time to move per cycle.  At high
    frequency it collapses to the straight line ``V/R(ic)``, and
    `test_memristor_loop_collapses_at_high_frequency` measures exactly
    that against the closed form.

    Biolek's window is what keeps ``x`` inside [0,1] *and* releases it
    again: ``stp(-i)`` is 1 when the current is negative, so the window
    that is ``1 - x^(2p)`` while the state is being driven up becomes
    ``1 - (x - 1)^(2p)`` when the current reverses.  Joglekar's earlier
    ``1 - (2x - 1)^(2p)`` vanishes at BOTH ends regardless of current
    direction, so a state that reaches a boundary is stuck there
    forever; that is the "boundary lock" problem and it is why the
    current-dependent form is the one worth having.

    **This model is the reason to want a state primitive the DSL does
    not have.**  Its state equation is ``dx/dt = g(x, ...)`` -- the
    window depends on the state being integrated -- and `idt` cannot
    express that, because an `idt` application cannot appear inside its
    own argument.  What is written instead is::

        xn = Node('xs'); bx = Branch(xn, minus, 'xdot')
        x  = idt(bx.V, ic)                 # the state, read back
        ...
        Contribution(bx.V, muv*ron/d**2 * i * w)   # the state equation

    an internal node whose potential *is* ``dx/dt``, driven by a
    potential contribution and integrated by `idt`.  Nothing else
    touches ``xn``, so KCL forces the branch current to zero and the
    node injects nothing into the circuit; it is a scratch signal
    (`ComparatorHdl` uses the same idiom for the latch).  The cost is
    one extra node, one extra branch current and a reader who has to be
    told why.  A `state(name, ic)` returning a symbol, plus a statement
    form for its equation, would be four lines and no explanation.

    **And it cannot use `var()` anywhere**, which is the second half of
    the same gap: the compiler collects `idt`/`idtmod`/`laplace_*`
    applications by walking the *statements*, never the let-chain, so an
    `idt` that appears only inside a `var()` definition never gets a
    state allocated and the model dies in the printer with
    ``Unsupported by _P: idt``.  `$limit` had exactly this bug and it
    was fixed (`hdl.py`, "Intermediates FIRST, and this is not
    tidiness"); the three state operators still have it.  Every
    intermediate of this model mentions ``x``, so the model is written
    flat -- which it can afford to be, and a MOSFET could not.

    The resistance reads a CLAMPED ``x``; the window reads the raw one.
    ``ron*x + roff*(1-x)`` passes through zero at ``x = roff/(roff-ron)``
    = 1.006 for the default card, six thousandths outside the physical
    range and a division by zero when a Newton iterate visits it.
    ``minc(maxc(x, 0), 1)`` is a hard clamp on an ATOM (the state symbol
    itself -- `differentiable-numerics`: never on a compound), it is
    exact everywhere the physics is defined, and outside [0,1] a
    constant resistance is the correct saturated answer anyway.  The
    window is deliberately left unclamped, because outside [0,1] it goes
    negative and pushes the state back -- checked on all four
    (sign of ``x``-excursion) x (sign of ``i``) combinations, it is
    restoring in every one.

    ``ic`` is spelled ``ic`` and not ``x0`` for a reason that is not
    taste: `Transient`'s ``uic`` seeding looks up an instance parameter
    **literally named ``ic``** before it will call the generated
    ``state_ic()``, so with any other spelling the state silently starts
    at zero while the DC pin (which reads the compiled ic function) is
    perfectly correct.  The symptom is a transient that begins from the
    wrong resistance and no error anywhere.
    """
    instparams = [Parameter(name='ron', desc='Resistance of the fully doped '
                            'film', unit='ohm', default=100.0),
                  Parameter(name='roff', desc='Resistance of the undoped '
                            'film', unit='ohm', default=16e3),
                  Parameter(name='muv', desc='Dopant mobility',
                            unit='m^2/(V*s)', default=1e-14),
                  Parameter(name='d', desc='Film thickness', unit='m',
                            default=10e-9),
                  Parameter(name='ic', desc='Initial doped fraction x(0)',
                            unit='', default=0.1),
                  Parameter(name='p', desc='Biolek window exponent',
                            unit='', default=1.0)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        xn = Node('xs')
        bx = Branch(xn, minus, 'xdot')
        ## The state.  `ic` pins it at dc (the LRM's "the initial
        ## condition shall force the DC solution") and seeds `uic`.
        x = idt(bx.V, ic)                                              # noqa
        xr = _minc(_maxc(x, 0.0), 1.0)
        i = b.V / (ron * xr + roff * (1 - xr))                         # noqa
        ## Biolek's `stp(-i)`, as a Piecewise on the SOLUTION.  Both
        ## arms are polynomials in a bounded quantity, so both are
        ## finite at every bias and the both-arms rule costs nothing
        ## here.  The window is discontinuous in `i` at `i = 0` and the
        ## drive `i*f` is not: the drive is C0 through the reversal,
        ## which is what the trajectory needs, and the jump lands in the
        ## Jacobian where Newton can see it rather than in the answer.
        stp = sympy.Piecewise((sympy.Integer(1), i < 0),
                              (sympy.Integer(0), True))
        ## `safe_abs(u)**(2*p)` rather than `u**(2*p)`: the exponent is a
        ## PARAMETER, so sympy cannot know it is even, and `u**(2*p)` is
        ## NaN for every negative `u` -- which is half the range of
        ## `x - stp`.  For an even power the absolute value is an exact
        ## identity, not an approximation, and `safe_abs` keeps the base
        ## strictly positive so the derivative is finite at `u = 0` too.
        w = 1 - _safe_abs(x - stp) ** (2 * p)                          # noqa
        return (Contribution(b.I, i),
                Contribution(bx.V, muv * ron / d ** 2 * i * w))        # noqa
