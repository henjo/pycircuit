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
``GummelPoonNpnHdl``  the SPICE Gummel-Poon bipolar transistor: one
``GummelPoonPnpHdl``  transport current divided by the base charge
                 ``qb``, which carries both Early voltages and both
                 high-injection knees; four base-current components;
                 base resistance modulated by ``qb``; both junctions'
                 depletion and diffusion charge with the ``xcjc``
                 split across the base resistance; the SPICE
                 temperature path; shot, flicker and thermal noise --
                 and the first model with TWO ``limit_pnj`` probes
                 sharing a terminal, which is the case the
                 write-back rule exists for
``MosLevel1Hdl``  the SPICE level-1 (Shichman-Hodges) MOSFET, and
``MosLevel1PmosHdl``  the first production user of ``limit_together``:
                 SPICE's own four limiter probes -- ``fetlim(vgs)``,
                 ``limvds(vds)``, ``pnjlim(vbs)``, ``pnjlim(vbd)`` --
                 which the per-probe form CANNOT COMPILE, because the
                 fourth finds both of its terminals claimed.  Also the
                 first model to use ``$param_given``, ``$limit`` and
                 `Collapse` together, and the one whose source/drain
                 swap is written without a ``sign``
``OpAmpHdl``     a Boyle-class opamp macromodel: input offsets and bias
                 currents, CMRR and PSRR, one dominant pole, slew-rate
                 limiting, rail clamping and output current limiting --
                 the last two from one ``tanh`` whose slope is
                 ``1/rout`` and whose asymptote is ``isc``.  Every gain
                 relation is an identity of Boyle's parameterisation,
                 which is what makes them assertable rather than
                 fittable
``GummelPoonNpnThermalHdl``  the self-heating bipolar -- `_gp_core`
                 with the thermal node bolted on, and thermal runaway
                 with its onset measured against ``rth*dP/dT = 1``
``EkvNmosHdl``   EKV 2.6's long-channel core: weak, moderate and
``EkvPmosHdl``   strong inversion from ONE interpolation function,
                 the four terminal charges with the Ward-Dutton
                 partition, and the first production user of
                 ``limit_fet``, ``limit_vds`` and ``softplus``.  An
                 independent formulation of what
                 `compact.PspMosLongChannel` solves exactly, which is
                 what makes the cross-check in
                 ``test_elements_hdl_library3.py`` worth having
``VcoHdl``       a VCO whose noise is injected into the FREQUENCY
                 through an ``idtmod`` phase, so jitter random-walks
``DividerHdl``   continuous divide-by-N, the second production user
                 of ``Cross``
``MixerHdl``     a multiplying mixer
``ChargePumpHdl``  UP/DN current pump for a PLL
``PhotodiodeHdl``  the SPICE diode with a photocurrent and an optical
``LedHdl``       node; solar-cell I-V and LED output linear above
                 threshold
``MesfetCurticeHdl``  Curtice 1980, Statz (``mes1.va``) and Angelov
``MesfetStatzHdl``    1992: tanh-shaped channels with pnjlim on the
``HemtAngelovHdl``    Schottky gates
``MosLevel3Hdl`` SPICE MOS level 3 (n and p), a transcription of
``MosLevel3PmosHdl``  ``mos3load.c`` with every knob; = level 1 at
                 gamma = 0.  All eleven use ``params_as = 'p'``
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
                                   laplace_nd, laplace_zp, TIME, TEMP,
                                   limit_pnj, limit_fet, limit_vds,
                                   limit_delta, limit_identity,
                                   limit_together,
                                   param_given)

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
from pycircuit.circuit.hdl import select as _select
from pycircuit.circuit.hdl import safe_pow as _safe_pow
from pycircuit.circuit.hdl import safe_sqrt as _safe_sqrt
from pycircuit.circuit.hdl import softplus as _softplus
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

       **The external network is in PARALLEL with this block's own
       ``rth``, not in series with it** (measured 2026-08-25,
       ``test_the_thermal_pin_takes_an_external_foster_ladder``), because
       the RC below and the dissipation both sit on the same
       ``(th, tha)`` branch.  So the obvious way to say "no internal
       thermal resistance, use my ladder" -- ``rth = 0`` -- does the
       opposite: it collapses the branch onto a zero-volt source and
       SHORTS the ladder out, which is the trap the last paragraph of
       this docstring describes for the shared-node case.  The only
       route today is a large ``rth``, and its value then leaks into the
       answer: ``rth = 1e6`` against a 5000 K/W ladder still carries
       0.5% of the heat.  The fix is small and belongs here rather than
       in the DSL -- an ``external=True`` that simply does not emit the
       RC contribution, chosen by a build-time Python ``if``, which is
       legal because it selects which expression to compile rather than
       which arm to evaluate.  Not done; recorded with its number so it
       can be measured.

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


def _spice_diode(p, a, c, T, junction=None):
    """The junction itself, shared by the isothermal and self-heating
    classes and by the optical devices.  Returns ``(statements,
    p_dissipated, i_junction)``.

    ``junction``, if given, is a callable ``junction(bd) -> expression``
    that is contributed as an EXTRA current across the junction branch
    ``bd`` (between the internal node behind ``rs`` and the cathode):
    the photodiode's photocurrent and shunt belong there, in parallel
    with the junction and behind the series resistance, and the branch
    is local to this function.  2026-08-26, fifth batch: a helper's
    internal node cannot be reached from the calling ``analog()``, so
    the helper has to accept the caller's contribution rather than the
    caller reaching in.  With ``junction=None`` nothing is added and
    the generated code is what it was.

    ``p`` is the calling class's parameter namespace (``params_as =
    'p'``; `hdl.ParamNamespace`).  Until 2026-08-26 this function took
    all nineteen parameters as explicit arguments, because a helper has
    its own ``__globals__`` and the bare names the compiler injects into
    ``analog()`` do not reach it; the namespace is roadmap S5's answer.

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
    tr = _var(T / p.tnom, 'tratio')
    ltr = _var(sympy.log(tr), 'ltratio')
    egT = _var(1.16 - 7.02e-4 * T ** 2 / (T + 1108.0), 'egT')
    egn = _var(1.16 - 7.02e-4 * p.tnom ** 2 / (p.tnom + 1108.0), 'egtnom')
    ## `expl`, not `exp`.  The argument is bounded ABOVE by
    ## `eg*q/(n*k*tnom)` (~43 for a silicon card at n=1), so on any
    ## sensible card the two are identical; `expl` costs nothing and
    ## keeps a card with n < 0.5 -- or a runaway thermal node -- finite.
    isT = _var(p.IS * _expl((tr - 1) * p.eg / (p.n * vtT)
                            + p.xti / p.n * ltr), 'isT')
    vjT = _var(p.vj * tr - 3 * vtT * ltr - egn * tr + egT, 'vjT')
    cjT = _var(p.cjo * (1 + p.m * (4e-4 * (T - p.tnom) - (vjT / p.vj - 1))),
               'cjT')

    ## -- static current ------------------------------------------------
    vd = _var(bd.V, 'vd')
    nvt = _var(p.n * vtT, 'nvt')
    ifwd = _var(isT * (_expl(vd / nvt) - 1), 'ifwd')
    ibrk = _var(p.ibv * _expl(-(vd + p.bv) / vtT), 'ibrk')
    idio = _var(p.area * (ifwd - ibrk), 'id')

    ## -- charge --------------------------------------------------------
    fcvj = _var(p.fc * vjT, 'fcvj')
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
    ufl = _var(0.5 * (1 - p.fc), 'ufloor')
    ubase = _var(ufl + _hypsmooth(ubar - ufl, 1e-9 * (1 - p.fc)), 'ubase')
    qlo = _var(cjT * vjT * (1 - ubase ** (1 - p.m)) / (1 - p.m), 'qdeplo')
    f1 = _var(vjT * (1 - (1 - p.fc) ** (1 - p.m)) / (1 - p.m), 'f1')
    f2 = _var((1 - p.fc) ** (1 + p.m), 'f2')
    f3 = _var(1 - p.fc * (1 + p.m), 'f3')
    qhi = _var(cjT * (f1 + (f3 * (vd - fcvj)
                            + p.m / (2 * vjT) * (vd ** 2 - fcvj ** 2)) / f2),
               'qdephi')
    qdep = _var(sympy.Piecewise((qlo, vd < fcvj), (qhi, True)), 'qdep')
    ## Diffusion charge from the FORWARD current only: the breakdown
    ## term is not minority-carrier injection and carries no transit
    ## time, and `tt*ibrk` at 1 kV reverse would be a large spurious
    ## capacitance.
    qtot = _var(p.area * (qdep + p.tt * ifwd), 'qtot')

    ## -- noise ---------------------------------------------------------
    ## `safe_abs`, not `Abs`: the PSD must be non-negative and the
    ## current is negative over most of the reverse range.  (`CY` is
    ## evaluated, never differentiated, so the regularisation here is
    ## about the value only.)  With `kf = 0` the flicker term is zero and
    ## still evaluated -- which is fine, because `|I|^af` is finite for
    ## every bias, floored below by `safe_abs`'s own eps.
    iabs = _var(_safe_abs(idio), 'iabs')
    noise = (Contribution(bd.I, _white_noise(2 * _QE * iabs)),
             Contribution(bd.I, _flicker_noise(p.kf * iabs ** p.af, 1)),
             ## Thermal noise of the series resistance.  It goes away
             ## with the branch when rs collapses, which is the correct
             ## behaviour and comes for free.
             Contribution(brs.I, _white_noise(4 * _KB * T * p.area / p.rs)))

    ## -- statements ----------------------------------------------------
    stmts = noise + (
        Contribution(brs.I, brs.V * p.area / p.rs),
        ## The Verilog-A optional-parasitic idiom, and the reason
        ## `rs = 0` (SPICE's default) is not a division by zero: the
        ## collapsed variant never compiles the `1/rs`.
        Collapse(brs, p.rs <= 0),
        Contribution(bd.I, idio),
        Contribution(bd.I, ddt(qtot)))
    if junction is not None:
        stmts = stmts + (Contribution(bd.I, junction(bd)),)
    ## Total static dissipation: the series chain carries one current, so
    ## the terminal-to-terminal voltage times it is the whole of it.  The
    ## charge term is deliberately absent -- storage does not dissipate.
    return stmts, _var(Branch(a, c).V * idio, 'pdiss'), idio
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
    params_as = 'p'
    instparams = _spice_diode_params()

    @staticmethod
    def analog(p, plus, minus):
        return _spice_diode(p, plus, minus, TEMP)[0]


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
    params_as = 'p'
    instparams = _spice_diode_params() + _thermal_params()

    @staticmethod
    def analog(p, plus, minus, th, tha):
        heat = SelfHeating(th, tha, p.rth, p.cth)
        stmts, pdiss, _idio = _spice_diode(p, plus, minus, heat.T)
        return stmts + heat.dissipate(pdiss)


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

    **It used to be unable to use `var()` anywhere**, and that is why it
    is written flat.  The compiler collected `idt`/`idtmod`/`laplace_*`
    applications by walking the *statements*, never the let-chain, so an
    `idt` appearing only inside a `var()` definition never got a state
    allocated and the model died in the printer with ``Unsupported by
    _P: idt``.  `$limit` had exactly the same bug and had already been
    fixed (`hdl.py`, "Intermediates FIRST, and this is not tidiness");
    the three state operators had not.

    **FIXED 2026-08-25** -- both state-collection loops now walk
    intermediates as well as statements, and this model reported the
    gap that prompted it.  The flat spelling below is therefore no
    longer forced, and rewriting it with a let-chain is a reasonable
    tidy-up; it is left flat because every intermediate here mentions
    ``x`` anyway, so the chain would buy nothing but churn on a model
    whose whole purpose is to be read.  A MOSFET could not have
    afforded the flat form, which is what made the bug worth fixing
    rather than documenting.

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


## ======================================================================
## The SPICE Gummel-Poon BJT -- roadmap item 3
## ======================================================================


def _pn_depletion_charge(v, cj, vjx, mx, fcx, tag):
    """SPICE's junction depletion charge, with the ``fc`` linearisation.

    The same construction the diode above uses, written once more --
    deliberately, and the difference is the point.  `_spice_diode`
    builds the floored power base BY HAND (three lines, a floor-placement
    argument and a smoothing-width scan) because `hdl.safe_pow` did not
    exist when it was written, and its comment says so and keeps the
    hand-rolled version as a record of the friction.  This is what the
    same block costs now::

        _safe_pow(1 - v/vj, 1 - m, lo=0.5*(1 - fc))

    one call, and the floor argument moves from a paragraph in the model
    to the kernel's own docstring.  The floor is at HALF the seam value
    ``1 - fc``, i.e. strictly inside the arm that is always discarded
    above the seam, so `safe_pow`'s hard `maxc` clamp is the right
    instrument: it is never reached by a selected evaluation.

    ``fc`` and ``m`` are clamped to SPICE's own limits (`fc < 1`,
    `m < 1`), because the model divides by ``1 - fc`` and by ``1 - m``.
    ngspice does this in `bjttemp.c`/`diotemp.c` with an outright
    assignment and a warning; here it is a `minc`, which is the same
    number for every card that is not already broken.
    """
    ## Parameter-only clamps: no bias reaches them, so their derivative
    ## with respect to a solution variable is exactly zero and the hard
    ## corner costs nothing.
    fcc = _var(_minc(fcx, 0.9999), 'fc' + tag)
    mc = _var(_minc(mx, 0.9), 'm' + tag)
    fcvj = _var(fcc * vjx, 'fcvj' + tag)
    ub = _var(_safe_pow(1.0 - v / vjx, 1.0 - mc, lo=0.5 * (1.0 - fcc)),
              'ub' + tag)
    qlo = _var(cj * vjx * (1.0 - ub) / (1.0 - mc), 'qlo' + tag)
    f1 = _var(vjx * (1.0 - (1.0 - fcc) ** (1.0 - mc)) / (1.0 - mc),
              'f1' + tag)
    f2 = _var((1.0 - fcc) ** (1.0 + mc), 'f2' + tag)
    f3 = _var(1.0 - fcc * (1.0 + mc), 'f3' + tag)
    qhi = _var(cj * (f1 + (f3 * (v - fcvj)
                           + mc / (2.0 * vjx) * (v * v - fcvj * fcvj)) / f2),
               'qhi' + tag)
    return _var(sympy.Piecewise((qlo, v < fcvj), (qhi, True)), 'qj' + tag)


def _spice_bjt_params():
    """The SPICE Gummel-Poon card, in SPICE's own names and defaults.

    ``is`` is a Python keyword, so the saturation current keeps the
    file's spelling ``IS``; everything else is lower-cased SPICE.

    Two of these names are traps worth stating rather than discovering.
    ``var`` is SPICE's REVERSE Early voltage and it collides by name
    with `hdl.var`, the let-chain binder -- which is precisely why this
    module imports that as ``_var`` (see the import block: the alias was
    kept for readability after the bug that forced it was fixed, and
    here is the case that would have needed it anyway).  Since
    2026-08-26 the BJT classes read their card through a parameter
    namespace (``p.var``), so the collision is confined to models that
    read parameters as bare names.  ``nc`` is the base-collector
    leakage emission coefficient, not a node count.

    SPICE's "0 means infinite" convention is kept for ``vaf``, ``var``,
    ``ikf`` and ``ikr``: a zero Early voltage or knee current has no
    physical reading, and every card in existence uses the convention.
    """
    return [
        Parameter(name='IS', desc='Transport saturation current at tnom',
                  unit='A', default=1e-16),
        Parameter(name='bf', desc='Ideal maximum forward beta', unit='',
                  default=100.0),
        Parameter(name='nf', desc='Forward current emission coefficient',
                  unit='', default=1.0),
        Parameter(name='vaf', desc='Forward Early voltage (0 = infinite)',
                  unit='V', default=0.0),
        Parameter(name='ikf', desc='Forward knee current (0 = infinite)',
                  unit='A', default=0.0),
        Parameter(name='ise', desc='Base-emitter leakage saturation '
                  'current', unit='A', default=0.0),
        Parameter(name='ne', desc='Base-emitter leakage emission '
                  'coefficient', unit='', default=1.5),
        Parameter(name='br', desc='Ideal maximum reverse beta', unit='',
                  default=1.0),
        Parameter(name='nr', desc='Reverse current emission coefficient',
                  unit='', default=1.0),
        Parameter(name='var', desc='Reverse Early voltage (0 = infinite)',
                  unit='V', default=0.0),
        Parameter(name='ikr', desc='Reverse knee current (0 = infinite)',
                  unit='A', default=0.0),
        Parameter(name='isc', desc='Base-collector leakage saturation '
                  'current', unit='A', default=0.0),
        Parameter(name='nc', desc='Base-collector leakage emission '
                  'coefficient', unit='', default=2.0),
        Parameter(name='rb', desc='Zero-bias base resistance', unit='ohm',
                  default=0.0),
        ## A NEGATIVE default as the "unset" sentinel, and it is a
        ## workaround for a live defect rather than a preference.
        ## SPICE's ``RBM`` defaults to ``RB``, which a default VALUE
        ## cannot express -- a card that deliberately says ``rbm = 0``
        ## wants a base resistance that collapses to nothing at high
        ## injection, and a card that says nothing wants a constant
        ## ``rb``.  `hdl.param_given` is exactly the instrument for
        ## that, and this model was written with it.
        ##
        ## THE REASON BELOW HAS EXPIRED; the rule has not.  This used
        ## to say `param_given` could not be used here because
        ## `$param_given` and `$limit` in one model broke the generated
        ## `limit()` at run time -- the parameter expressions were
        ## lambdified over `paramsyms + [TEMP]` while `limit()` called
        ## them with the givenness flags appended, one argument too
        ## many.  **That was fixed on 2026-08-25**
        ## (`test_limit_with_param_given.py`), and `MosLevel1Hdl` below
        ## uses both together for exactly this kind of rule.
        ##
        ## The sentinel is KEPT anyway, and deliberately: `rbm = -1.0`
        ## is the shipped default and an explicit `rbm = -1` card would
        ## change meaning if this became `param_given`.  That is a
        ## choice now rather than a workaround, and the cost of the
        ## choice is the sentence below -- `rbm = -1` and "rbm not
        ## given" are the same thing here and are two different things
        ## in SPICE.
        ##
        ## A negative resistance is not a physical card, so it is a
        ## sound sentinel -- but it is a sentinel, and the difference
        ## matters: `rbm = -1` and "rbm not given" are the same thing
        ## here and are two different things in SPICE.
        Parameter(name='rbm', desc='Minimum base resistance at high '
                  'current (negative = follow rb)', unit='ohm',
                  default=-1.0),
        Parameter(name='re', desc='Emitter resistance', unit='ohm',
                  default=0.0),
        Parameter(name='rc', desc='Collector resistance', unit='ohm',
                  default=0.0),
        Parameter(name='cje', desc='Zero-bias base-emitter depletion '
                  'capacitance', unit='F', default=0.0),
        Parameter(name='vje', desc='Base-emitter built-in potential',
                  unit='V', default=0.75),
        Parameter(name='mje', desc='Base-emitter grading coefficient',
                  unit='', default=0.33),
        Parameter(name='tf', desc='Ideal forward transit time', unit='s',
                  default=0.0),
        Parameter(name='xtf', desc='Bias coefficient of tf', unit='',
                  default=0.0),
        Parameter(name='vtf', desc='Vbc dependence of tf (0 = none)',
                  unit='V', default=0.0),
        Parameter(name='itf', desc='High-current parameter of tf',
                  unit='A', default=0.0),
        Parameter(name='cjc', desc='Zero-bias base-collector depletion '
                  'capacitance', unit='F', default=0.0),
        Parameter(name='vjc', desc='Base-collector built-in potential',
                  unit='V', default=0.75),
        Parameter(name='mjc', desc='Base-collector grading coefficient',
                  unit='', default=0.33),
        Parameter(name='xcjc', desc='Fraction of cjc connected to the '
                  'INTERNAL base node', unit='', default=1.0),
        Parameter(name='tr', desc='Ideal reverse transit time', unit='s',
                  default=0.0),
        Parameter(name='fc', desc='Forward-bias depletion coefficient',
                  unit='', default=0.5),
        Parameter(name='xtb', desc='Forward and reverse beta temperature '
                  'exponent', unit='', default=0.0),
        Parameter(name='eg', desc='Energy gap', unit='eV', default=1.11),
        Parameter(name='xti', desc='Saturation-current temperature '
                  'exponent', unit='', default=3.0),
        Parameter(name='kf', desc='Flicker-noise coefficient', unit='',
                  default=0.0),
        Parameter(name='af', desc='Flicker-noise exponent', unit='',
                  default=1.0),
        Parameter(name='area', desc='Area scaling factor', unit='',
                  default=1.0),
        Parameter(name='tnom', desc='Parameter measurement temperature',
                  unit='K', default=300.15),
    ]



def _gp_core(p, T, npn, c, b, e):
    """The Gummel-Poon body itself, shared by the isothermal and the
    self-heating classes.  Returns ``(statements, (ict, ibc, ibe))``.

    ``p`` is the calling class's parameter namespace
    (`hdl.ParamNamespace`), and every one of the card's forty-two
    parameters is read as ``p.<name>``.  Until 2026-08-26 this function
    read them as bare names it did not define, and was callable only
    through a hand-applied ``types.FunctionType`` rebind of its globals
    onto the calling ``analog()``'s injected copy (``_with_params``):
    the compiler binds a class's parameter symbols into ``analog()``'s
    OWN globals (`hdl._analog_function`), so a helper shared by two
    classes could not see them, and `_spice_diode` paid for the same
    thing with a nineteen-argument signature.  Roadmap S5's namespace
    is the fix, and it is bit-identical: ``p.bf`` IS the symbol the
    bare style bound (``test_hdl_params.py`` pins the three BJT classes
    to hashes recorded before the change).

    The two ``limit_pnj`` declarations read the DEVICE'S OWN ``isT`` and
    ``vtT`` -- the heated ones, for the self-heating variant.  Until
    2026-08-25 that was refused ("a limiter's parameters are evaluated
    BEFORE the device is, so they may not reach the solution"), and this
    function took a ``tlim`` argument so the thermal class could hand
    the limiter an AMBIENT copy of the saturation current -- five
    duplicated lines, and a critical voltage placed against a junction
    up to 100 K hotter than it thought.  The refusal was right about
    the order and wrong about the conclusion: a limiter is handed the
    last accepted iterate precisely so it can measure against it, and
    a parameter evaluated THERE is as well-defined as ``vold`` is.
    `hdl._limit_par_fn` now does exactly that, and ``tlim`` is gone.

    The equations, the references and what is deliberately absent are in
    `_gummel_poon`.
    """
    ## Three internal nodes, one per parasitic resistance, each
    ## collapsed away when its resistance is zero -- which is
    ## SPICE's default for all three, so a bare instance has no
    ## internal nodes at all and the `1/rb` is never compiled.
    bi, ci, ei = Node('bi'), Node('ci'), Node('ei')
    brb, brc, bre = Branch(b, bi), Branch(c, ci), Branch(e, ei)
    ## POLARITY LIVES IN THE BRANCH TERMINAL ORDER, and nowhere
    ## else.  A p-n-p's branches are declared reversed, so `bbe.V`
    ## is the n-p-n-convention forward bias for both devices and
    ## every expression below is written once; a `Contribution` on
    ## the reversed branch then delivers the current in the
    ## reversed direction too, which is exactly the polarity flip
    ## SPICE applies term by term on the way out.
    ##
    ## Doing it this way rather than with a `type` factor is not
    ## only tidier: `limit_pnj` limits a BRANCH POTENTIAL in the
    ## solution vector, so a p-n-p whose junction voltage is
    ## `-V(bi,ei)` would have the limiter bound the wrong sign --
    ## compressing reverse excursions and letting forward ones
    ## through.  Reversing the branch puts the limiter on the
    ## quantity that actually has the exponential in it.
    ##
    ## The two junctions SHARE THE INTERNAL BASE as one terminal.
    ## For the n-p-n it is the PLUS of both, which is the case
    ## `limit_spec`'s write-back rule exists for: the first probe
    ## moves `bi`, so the second must move its own minus (`ci`) or
    ## it would undo the first.  For the p-n-p the base is the
    ## minus of both, the plus terminals differ, and each probe
    ## moves its own -- the same rule, taking its other branch.
    ##
    ## The transport current is its own branch from internal
    ## collector to internal emitter -- ONE current source, not two.
    ##
    ## The external fraction of the collector capacitance
    ## (`1 - xcjc`) hangs off the EXTERNAL base, behind rb.  That is
    ## the whole point of the parameter: the collector-base overlap
    ## capacitance is distributed along the base resistance.
    if npn > 0:
        bbe, bbc, bct = Branch(bi, ei), Branch(bi, ci), Branch(ci, ei)
        bbx = Branch(b, ci)
    else:
        bbe, bbc, bct = Branch(ei, bi), Branch(ci, bi), Branch(ei, ci)
        bbx = Branch(ci, b)

    ## -- temperature ------------------------------------------------
    vtT = _var(_vt(T), 'vtT')
    trat = _var(T / p.tnom, 'trat')
    ltr = _var(sympy.log(trat), 'ltrat')
    ## SPICE's `factlog`, `bjttemp.c`: one exponent shared by the
    ## transport current and (scaled by 1/NE, 1/NC) the two leakage
    ## currents.  `expl` rather than `exp` for the same reason the
    ## diode above uses it -- the argument is bounded above on any
    ## sensible card, so the two are the same function there, and a
    ## card with a small emission coefficient stays finite.
    factlog = _var((trat - 1.0) * p.eg / vtT + p.xti * ltr, 'factlog')
    isT = _var(p.area * p.IS * _expl(factlog), 'isT')
    ## `XTB` moves both betas and, inversely, both leakage currents
    ## -- SPICE divides the leakage by the same factor it multiplies
    ## beta by, which is what keeps the low-injection Gummel plot's
    ## SHAPE roughly fixed with temperature.
    bfac = _var(_expl(p.xtb * ltr), 'bfac')
    bfT = _var(p.bf * bfac, 'bfT')
    brT = _var(p.br * bfac, 'brT')
    iseT = _var(p.area * p.ise * _expl(factlog / p.ne) / bfac, 'iseT')
    iscT = _var(p.area * p.isc * _expl(factlog / p.nc) / bfac, 'iscT')
    ## Junction potentials and capacitances: the diode's path, per
    ## junction.  `egn` is the gap at tnom, `egT` at T.
    egT = _var(1.16 - 7.02e-4 * T ** 2 / (T + 1108.0), 'egT')
    egn = _var(1.16 - 7.02e-4 * p.tnom ** 2 / (p.tnom + 1108.0),
               'egtnom')
    vjeT = _var(p.vje * trat - 3.0 * vtT * ltr
                - egn * trat + egT, 'vjeT')
    vjcT = _var(p.vjc * trat - 3.0 * vtT * ltr
                - egn * trat + egT, 'vjcT')
    cjeT = _var(p.area * p.cje * (1.0 + p.mje * (4e-4 * (T - p.tnom)
                                           - (vjeT / p.vje - 1.0))),
                'cjeT')
    cjcT = _var(p.area * p.cjc * (1.0 + p.mjc * (4e-4 * (T - p.tnom)
                                           - (vjcT / p.vjc - 1.0))),
                'cjcT')

    ## -- the two junction voltages, LIMITED ------------------------
    ## This used to spell the temperature-scaled saturation current
    ## TWICE -- once as `isT` in the chain, once bare here -- because
    ## `limit_pnj`'s parameters were lambdified over the parameters
    ## and TEMP alone, so an intermediate symbol was not in scope.
    ## Roadmap §12.4 fixed that: the compiler now resolves a `var()`
    ## symbol against the chain (reading it, not inlining it), so the
    ## natural spelling works and the duplicate is gone.
    ##
    ## ngspice computes its critical voltage from the UNSCALED
    ## saturation current (`bjttemp.c` sets `tVcrit` before
    ## `bjtload.c` multiplies by area), and uses the same one for
    ## both junctions.  Passing the area-scaled current instead
    ## moves `vc` by `VT*ln(area)` -- 26 mV per decade of area, on a
    ## quantity that only decides where compression BEGINS.  The
    ## scaled one is the device's actual saturation current, so it
    ## is the one used here; the difference is recorded because it
    ## is a real, if small, divergence from the reference.
    ## The heated `isT`/`vtT` themselves: a parameter that reaches the
    ## solution is evaluated at the LAST ACCEPTED iterate (SPICE's `von`
    ## semantics), so the limiter sees the junction at the temperature
    ## it actually has.
    vbe = _var(limit_pnj(bbe.V, isT, p.nf * vtT), 'vbe')
    vbc = _var(limit_pnj(bbc.V, isT, p.nr * vtT), 'vbc')

    ## -- transport and base currents -------------------------------
    ifwd = _var(isT * (_expl(vbe / (p.nf * vtT)) - 1.0), 'ifwd')
    irev = _var(isT * (_expl(vbc / (p.nr * vtT)) - 1.0), 'irev')

    ## SPICE's "0 means infinite" for the Early voltages and the
    ## knee currents.  The reciprocal is taken inside the guarded
    ## arm AND floored, because both arms of a `Piecewise` are
    ## evaluated: `1/vaf` at `vaf = 0` is `inf`, and an `inf` in a
    ## discarded arm is only harmless until something multiplies it
    ## by the selected arm's zero derivative.
    ##
    ## `hdl.select` cannot take this one over, and the reason is worth
    ## the four lines: two of the four calls pass a PRODUCT
    ## (`area*ikf`), so the condition is `area*ikf <= 0`, which
    ## constrains neither factor on its own.  Clamping the product
    ## itself would need a substitution for the compound -- and sympy
    ## does not keep the compound: `1.0/(area*ikf)` is stored as
    ## `1.0 * area**-1 * ikf**-1`, in which `area*ikf` is not a node,
    ## so the substitution would silently do nothing.  Binding the
    ## product with `var()` first would give `select` an atom to work
    ## on; it also changes the emitted chain, so it is not done here.
    def _recip(p):
        return sympy.Piecewise((sympy.Float(0.0), p <= 0.0),
                               (1.0 / _maxc(p, 1e-30), True))
    rvaf = _var(_recip(p.vaf), 'rvaf')
    rvar = _var(_recip(p.var), 'rvar')
    rikf = _var(_recip(p.area * p.ikf), 'rikf')
    rikr = _var(_recip(p.area * p.ikr), 'rikr')

    ## THE BASE CHARGE.  `q1` is the Early (base-width modulation)
    ## factor and `q2` the high-injection one; `qb` is what every
    ## transport current is divided by.
    ##
    ## The floor on `q1`'s denominator is a domain guard, not a fit:
    ## it fires only when a junction is forward-biased to within 1%
    ## of an Early voltage -- 99 V on a 100 V card -- which no
    ## physical bias reaches and a wild Newton iterate does.
    q1d = _var(_maxc(1.0 - vbc * rvaf - vbe * rvar, 1e-2), 'q1d')
    q1 = _var(1.0 / q1d, 'q1')
    q2 = _var(ifwd * rikf + irev * rikr, 'q2')
    ## `1 + 4*q2` is >= 1 - 4*IS*(1/IKF + 1/IKR) by construction, so
    ## the floor is unreachable in the selected region and exists to
    ## keep the square root's argument in its domain for a Newton
    ## iterate that has driven both junctions hard into reverse.
    qb = _var(q1 * 0.5 * (1.0 + sympy.sqrt(_maxc(1.0 + 4.0 * q2,
                                                 1e-8))), 'qb')

    ict = _var((ifwd - irev) / qb, 'ict')
    ibe = _var(ifwd / bfT + iseT * (_expl(vbe / (p.ne * vtT))
                                    - 1.0), 'ibe')
    ibc = _var(irev / brT + iscT * (_expl(vbc / (p.nc * vtT))
                                    - 1.0), 'ibc')

    ## -- charge ----------------------------------------------------
    ## The transit-time modulation, SPICE's
    ## `tf_eff = TF*(1 + XTF*(IF/(IF + ITF))^2 * exp(VBC/(1.44*VTF)))`.
    ##
    ## `safe_div` rather than a division: with `itf = 0` -- SPICE's
    ## default, and the case where SPICE skips the factor entirely
    ## -- the ratio is `IF/IF`, and `safe_div` returns 1 for any
    ## current well above its epsilon and 0 at exactly zero
    ## current, where the charge it multiplies is zero anyway.  That
    ## reproduces SPICE's branch without a branch.
    tfr = _var(_safe_div(ifwd, ifwd + p.area * p.itf), 'tfr')
    ## `vtf = 0` means "no Vbc dependence".  The discarded arm is
    ## floored at 1 mV so that it cannot overflow while it is being
    ## discarded -- a card with `0 < vtf < 1e-3` would be clamped,
    ## and no such card exists.  That last sentence is the `margin`
    ## trade `select` makes the author state out loud: with
    ## `margin=0` the clamp would be `maxc(vtf, 0)` and the discarded
    ## arm would divide by zero.
    tfx = _var(_select(
        (sympy.Float(1.0), p.vtf <= 0.0),
        (_expl(vbc / (1.44 * p.vtf)), True), margin=1e-3), 'tfx')
    tff = _var(p.tf * (1.0 + p.xtf * tfr * tfr * tfx), 'tff')

    qbe = _var(tff * ifwd
               + _pn_depletion_charge(vbe, cjeT, vjeT, p.mje,
                                      p.fc, 'e'), 'qbe')
    qbc = _var(p.tr * irev
               + p.xcjc * _pn_depletion_charge(vbc, cjcT, vjcT,
                                             p.mjc, p.fc, 'c'),
               'qbc')
    ## The external fraction, on its own branch and therefore its
    ## own voltage: `V(b, ci)`, not `V(bi, ci)`.  With `rb = 0` the
    ## two are the same node pair and this is arithmetically the
    ## rest of `cjc`; with `rb > 0` it is the distributed part.
    qbx = _var((1.0 - p.xcjc)
               * _pn_depletion_charge(bbx.V, cjcT, vjcT, p.mjc,
                                      p.fc, 'x'), 'qbx')

    ## -- base resistance -------------------------------------------
    ## `rbm < 0` means "unset", i.e. `rbm = rb` (SPICE).  See the
    ## parameter's own comment for why this is a sentinel and not
    ## `param_given`, which is the right instrument and is currently
    ## incompatible with `$limit`.
    rbmx = _var(sympy.Piecewise(
        (p.rb, p.rbm < 0.0), (p.rbm, True)), 'rbmx')
    ## Base-width modulation: the spreading resistance falls as the
    ## base charge grows.  `qb >= 0.5*q1 > 0` structurally, so this
    ## is a division by something bounded away from zero.
    rbb = _var((rbmx + (p.rb - rbmx) / qb) / p.area, 'rbb')

    ## -- statements ------------------------------------------------
    ## Every expression above is in the n-p-n convention; the
    ## p-n-p's polarity is already carried by its reversed branch
    ## declarations, so nothing here is multiplied by a sign.
    stmts = (
        Contribution(bct.I, ict),
        Contribution(bbe.I, ibe + ddt(qbe)),
        Contribution(bbc.I, ibc + ddt(qbc)),
        Contribution(bbx.I, ddt(qbx)),
        Contribution(brb.I, brb.V / rbb),
        Collapse(brb, p.rb <= 0.0),
        Contribution(brc.I, brc.V * p.area / p.rc),
        Collapse(brc, p.rc <= 0.0),
        Contribution(bre.I, bre.V * p.area / p.re),
        Collapse(bre, p.re <= 0.0),
    )

    ## -- noise -----------------------------------------------------
    ## Shot noise is `2*q*|I|` on each of the two independent
    ## carrier streams -- the base current and the collector
    ## current -- and they are independent, which is why they are
    ## two uncorrelated contributions on two different branches and
    ## not one on the emitter.  `safe_abs` because a PSD may not be
    ## negative and both currents are negative over part of the
    ## range.
    icol = _var(_safe_abs(ict - ibc), 'icabs')
    ibas = _var(_safe_abs(ibe + ibc), 'ibabs')
    noise = (
        Contribution(bct.I, _white_noise(2.0 * _QE * icol)),
        Contribution(bbe.I, _white_noise(2.0 * _QE * ibas)),
        Contribution(bbe.I, _flicker_noise(p.kf * ibas ** p.af, 1)),
        ## Thermal noise of the three parasitics.  Each goes away
        ## with its branch when the resistance collapses, which is
        ## the correct behaviour and comes for free.  The base one
        ## is BIAS-DEPENDENT, because `rbb` is.
        Contribution(brb.I, _white_noise(4.0 * _KB * T / rbb)),
        Contribution(brc.I,
                     _white_noise(4.0 * _KB * T * p.area / p.rc)),
        Contribution(bre.I,
                     _white_noise(4.0 * _KB * T * p.area / p.re)),
    )
    ## The three static currents are returned so that a self-heating
    ## variant can build the dissipation from EXTERNAL terminal
    ## voltages.  Doing it that way rather than summing `V*I` over the
    ## internal branches is not a style choice: `brc.V**2*area/rc`
    ## carries a `1/rc`, and `rc = 0` is the default, so the parasitic
    ## terms have to be reached through voltages that already contain
    ## their IR drops.  See `GummelPoonNpnThermalHdl`.
    return stmts + noise, (ict, ibc, ibe)
def _gummel_poon(npn):
    """Build the ``analog()`` body of a Gummel-Poon BJT of one polarity.

    ``npn`` is +1 for an n-p-n and -1 for a p-n-p, closed over at class
    creation exactly as `compact._psp_mos_analog` closes over the
    channel type: the polarity cannot vary with bias, so branching on it
    at run time would double the compiled expression to no purpose.
    Every internal voltage below is ``npn*(terminal - terminal)``, so the
    body is written once for the n-p-n convention, and the polarity is
    restored on the contributions.

    Equations: Massobrio & Antognetti, *Semiconductor Device Modeling
    with SPICE*, 2nd ed., ch. 2 -- the integral charge-control model of
    Gummel & Poon (BSTJ 49, 827, 1970) as SPICE2 reduced it, which is
    what ngspice's `bjt` implements.  The pieces, in the order they are
    written:

    * forward and reverse transport currents ``IF``/``IR``;
    * the normalised base charge ``qb = q1/2*(1 + sqrt(1 + 4*q2))``,
      with ``q1`` carrying the two Early voltages and ``q2`` the two
      high-injection knees;
    * the transport current ``ICT = (IF - IR)/qb`` -- ONE current, which
      is what makes this a charge-control model rather than two diodes
      and a pair of alphas;
    * four base-current components: ideal forward ``IF/BF``, ideal
      reverse ``IR/BR``, and the two non-ideal recombination terms
      ``ISE``/``NE`` and ``ISC``/``NC``, which are what bend the Gummel
      plot's base current at low injection;
    * base, emitter and collector resistances on internal nodes, the
      base one MODULATED by ``qb`` (base-width modulation makes the
      base spreading resistance fall as the device turns on);
    * depletion charge on both junctions with SPICE's ``fc``
      linearisation, diffusion charge ``tf_eff*IF`` and ``tr*IR``, and
      the ``xcjc`` split of the collector capacitance between the
      internal and the external base node;
    * the SPICE temperature path for ``IS``, ``BF``/``BR``, the leakage
      currents, both junction potentials and both capacitances;
    * shot noise on the base and collector currents, flicker noise on
      the base current, thermal noise on all three resistances.

    **Not implemented, named rather than left to be discovered:** the
    substrate junction (``cjs``/``vjs``/``mjs``), which needs a fourth
    pin that would float on every card that does not give it -- use a
    `DiodeSpiceHdl` between collector and substrate; ``irb``, the
    current where base resistance falls to half, whose SPICE form is
    ``3*(rb - rbm)*(tan(z) - z)/(z*tan(z)^2)`` with a removable
    singularity at ``z = 0`` that no both-arms-safe form expresses
    cheaply (the ``qb`` modulation SPICE uses when ``irb`` is unset IS
    implemented, and is what nearly every card relies on); and ``ptf``,
    excess phase, which is a delay and therefore wants ``absdelay``.
    """
    def analog(p, c, b, e):
        return _gp_core(p, TEMP, npn, c, b, e)[0]
    return analog



class GummelPoonNpnHdl(Behavioural):
    """The SPICE Gummel-Poon n-p-n bipolar transistor.

    Terminals ``(c, b, e)``.  Parameters and defaults are SPICE's, so a
    bare ``GummelPoonNpnHdl('c', 'b', 'e')`` is an ideal transport
    device with ``bf = 100``, no parasitic resistance, no charge, no
    leakage and no high-injection roll-off, and every block is switched
    on by giving its parameter.

    What makes this the charge-control model rather than two diodes:
    there is ONE transport current ``(IF - IR)/qb`` between collector
    and emitter, and the base charge ``qb`` in its denominator carries
    both Early voltages and both high-injection knees at once.  Two
    consequences are worth stating because they are what the model is
    tested on:

    * **beta is not a parameter of the current**, it is a measured
      ratio.  ``bf`` sets the ideal base current ``IF/BF``; the measured
      ``Ic/Ib`` is ``BF/qb``, so it falls to ``BF/2`` exactly where
      ``qb = 2``, which is ``IF = 2*IKF`` -- i.e. at a COLLECTOR current
      of exactly ``IKF``.  That identity is what ``ikf`` means and it is
      what the Gummel plot measures;
    * **the Early effect is an output conductance**, ``go = Ic/(VAF +
      VCE - VBE)`` in forward active, which is a statement about the
      slope of the measured ``Ic(Vce)`` and not about a parameter.

    ``limit_pnj`` is declared on both junctions.  This is the first
    model in the tree with TWO limited probes sharing a terminal (the
    internal base).  Which node each probe moves is decided at RUN TIME
    -- the terminal that drifted further from the last accepted point,
    the shared base going to whichever junction applies the larger
    correction (roadmap 12.1) -- so neither probe can undo the other.
    (This used to describe a compile-time rule, "the base-emitter probe
    moves the base, so the base-collector probe moves the collector";
    that rule was replaced on 2026-08-25.)  ``explain()`` prints both.

    See `GummelPoonPnpHdl` for the p-n-p, which is this model with its
    branches declared the other way round and nothing else changed.
    """
    params_as = 'p'
    instparams = _spice_bjt_params()

    analog = staticmethod(_gummel_poon(+1))


class GummelPoonPnpHdl(Behavioural):
    """The SPICE Gummel-Poon p-n-p bipolar transistor.

    `GummelPoonNpnHdl` with every branch declared with its terminals
    exchanged, which is the whole of the difference: the equations are
    written once in the n-p-n convention and the polarity lives in the
    branch orientation (see `_gummel_poon`).  Parameters are the same
    card with the same signs -- a p-n-p's ``IS``, ``cje`` and ``vaf``
    are positive magnitudes exactly as an n-p-n's are, and the device
    conducts for ``V(e) > V(b) > V(c)``.

    Redefining ``analog`` here is REQUIRED and not decorative: the DSL
    keys compilation on the name appearing in the class body, so a
    subclass that merely inherited it would silently reuse the n-p-n
    expression and be wrong with nothing to report it.
    """
    params_as = 'p'
    instparams = _spice_bjt_params()

    analog = staticmethod(_gummel_poon(-1))


class GummelPoonNpnThermalHdl(Behavioural):
    """`GummelPoonNpnHdl` with the thermal node attached -- the
    self-heating bipolar transistor.  Terminals ``(c, b, e, th, tha)``.

    The diff against `GummelPoonNpnHdl` is the four lines below: build
    the `SelfHeating` block, pass its ``T`` in place of ``TEMP``, and
    append its statements.  ``rth`` and ``cth`` join the card; with
    ``rth`` at its default of 0 the thermal branch collapses to a
    zero-volt source and this element is `GummelPoonNpnHdl` again, to the
    last digit (``test_thermal_bjt_at_rth_zero_is_the_isothermal_bjt``).

    **Why the bipolar and not another diode.**  Thermal runaway is the
    classic bipolar failure mode and it is a genuine positive feedback
    loop, not a monotone drift: at a fixed base-emitter voltage the
    collector current rises with temperature (``IS`` goes up faster than
    ``VT`` widens the exponential), which raises the dissipation, which
    raises the temperature.  The loop closes, and the model has a DC
    solution only while the loop gain

        L = rth * dP/dT|_{Vbe, Vce}

    is below 1.  Above it there is no low-temperature solution at all.
    Because the thermal node is a real unknown with an exact Jacobian
    rather than an outer iteration, that onset is a property of the
    Newton system and can be measured: ``test_thermal_runaway_onset_is_
    where_the_loop_gain_reaches_one`` brackets the critical ``rth`` from
    the compiled model and compares it against ``L = 1`` computed from an
    independent transcription, and ``test_thermal_runaway_is_a_transient
    _that_does_not_settle`` shows the same device settling below the
    critical resistance and running away above it.

    **The dissipation is the STATIC dissipation**, built from the two
    external terminal voltages and the three static currents::

        P = V(c,e)*(ICT - IBC) + V(b,e)*(IBE + IBC)

    Two things follow, and both are deliberate.  Charge storage does not
    dissipate, so the displacement currents are excluded -- which means
    that during a fast transient the parasitic resistances' dissipation
    is understated, since the external voltages carry the IR drop of the
    total current while the currents here are the static part.  And the
    parasitic dissipation is included *without* any ``1/rc`` appearing in
    the expression, which matters because ``rc = 0`` is SPICE's default:
    summing ``V*I`` over the internal branches instead would put a
    division by the collapsed resistance into the thermal contribution,
    where no `Collapse` can remove it.

    A p-n-p version is `_gp_core(p, heat.T, -1, ...)` and four more lines;
    it is not shipped because nothing in the tree has asked for one and
    every class costs import time.
    """
    params_as = 'p'
    instparams = _spice_bjt_params() + _thermal_params()

    @staticmethod
    def analog(p, c, b, e, th, tha):
        heat = SelfHeating(th, tha, p.rth, p.cth)
        stmts, (ict, ibc, ibe) = _gp_core(p, heat.T, +1, c, b, e)
        pdiss = _var(Branch(c, e).V * (ict - ibc)
                     + Branch(b, e).V * (ibe + ibc), 'pdiss')
        return stmts + heat.dissipate(pdiss)



## ======================================================================
## EKV 2.6 -- roadmap item 9
## ======================================================================


def _ekv_params():
    """The EKV long-channel card.

    Names are the EKV 2.6 documentation's, lower-cased.  ``cox`` and
    ``kp`` are both here and are NOT redundant: ``kp = mu*cox`` sets the
    current and ``cox`` alone sets the charge, and a card may give one
    without the other.  ``cox = 0`` is the default and switches the
    charge model off entirely -- the charges are a PRODUCT with ``cox``,
    never a quotient, so zero really is zero and not an infinity waiting
    to be multiplied by one.
    """
    return [
        Parameter(name='vto', desc='Long-channel threshold voltage',
                  unit='V', default=0.5),
        Parameter(name='gamma', desc='Body effect factor',
                  unit='V^0.5', default=0.7),
        Parameter(name='phi', desc='Bulk Fermi potential (2*phi_F)',
                  unit='V', default=0.7),
        Parameter(name='kp', desc='Transconductance parameter mu*Cox',
                  unit='A/V^2', default=50e-6),
        Parameter(name='cox', desc='Gate oxide capacitance per area '
                  '(0 = no charge model)', unit='F/m^2', default=0.0),
        Parameter(name='theta', desc='Mobility reduction coefficient',
                  unit='1/V', default=0.0),
        Parameter(name='w', desc='Channel width', unit='m', default=1e-6),
        Parameter(name='l', desc='Channel length', unit='m', default=1e-6),
        Parameter(name='tcv', desc='Threshold voltage temperature '
                  'coefficient', unit='V/K', default=0.0),
        Parameter(name='bex', desc='Mobility temperature exponent',
                  unit='', default=0.0),
        Parameter(name='kf', desc='Flicker-noise coefficient', unit='',
                  default=0.0),
        Parameter(name='af', desc='Flicker-noise exponent', unit='',
                  default=1.0),
        Parameter(name='ef', desc='Flicker-noise frequency exponent',
                  unit='', default=1.0),
        Parameter(name='tnom', desc='Parameter measurement temperature',
                  unit='K', default=300.15),
    ]


def _ekv_analog(T, nmos):
    """Build the ``analog()`` body of an EKV long-channel MOSFET.

    Reference: C. Enz, F. Krummenacher and E. Vittoz, "An analytical MOS
    transistor model valid in all regions of operation and dedicated to
    low-voltage and low-current applications", *Analog Integrated
    Circuits and Signal Processing* **8**, 83-114 (1995), and the EKV
    v2.6 model documentation (Bucher, Lallement, Enz, Theodoloz and
    Krummenacher, EPFL, 1998).  The charge expressions are the
    charge-based formulation of Enz and Vittoz, *Charge-Based MOS
    Transistor Modeling* (Wiley, 2006), ch. 4.

    **The interpolation function is the whole point.**  One expression

        i_f = ln(1 + exp((VP - VS)/(2*UT)))^2

    covers weak, moderate and strong inversion with no regional
    equations and nothing pasted between them.  Its two limits are the
    two textbook laws and they are what the tests assert:

    * ``x << 0``: ``ln(1 + e^x) -> e^x``, so ``i_f -> e^{(VP-VS)/UT}``
      and the drain current is exponential in the gate voltage with the
      slope factor ``n`` -- ``ln(10)*n*UT`` per decade, the subthreshold
      swing;
    * ``x >> 0``: ``ln(1 + e^x) -> x``, so
      ``ID -> (n*beta/2)*[(VP - VS)^2 - (VP - VD)^2]``, the square law,
      and ``(n*beta/2)*(VP - VS)^2`` in saturation.

    The join is not a smoothing function: it is one analytic function
    whose asymptotes are those two, which is what "all-region" means and
    what distinguishes this family from a threshold-voltage model with
    a subthreshold patch.

    ``ln(1 + exp(x))`` is `hdl.softplus`, which is total: the literal
    form overflows at ``x = 710`` and the value it was about to return
    was 710.  This is that function's first production user.

    **Implemented:** the pinch-off voltage with its body effect, the
    bias-dependent slope factor, the all-region interpolation, simple
    mobility reduction (``theta``), the four terminal charges with the
    Ward-Dutton partition, the EKV temperature path for ``vto``, ``kp``
    and ``phi``, channel thermal noise from the Klaassen-Prins integral
    and flicker noise.

    **Not implemented, named rather than left to be discovered:**
    velocity saturation (``ucrit``), channel-length modulation
    (``lambda``), short- and narrow-channel charge sharing
    (``leta``/``weta``/``dl``/``dw``), reverse short-channel effect
    (``q0``/``lk``), impact ionisation (``iba``/``ibb``/``ibn``), gate
    and junction leakage, and the source/drain sheet resistance.  What
    is here is the long-channel core, which is the part that is worth
    cross-checking against a surface-potential model: everything in that
    list is a correction bolted on top of it, and including them would
    only make the comparison harder to read.
    """
    def analog(d, g, s, b):
        ## POLARITY IN THE BRANCH ORDER, exactly as `_gummel_poon` does
        ## it and for the same reason: `limit_fet` and `limit_vds` bound
        ## a branch potential in the SOLUTION vector, so a p-channel
        ## whose gate drive is `-V(g,s)` would have `fetlim` compress
        ## the wrong direction.  Declaring the branches reversed puts
        ## every probe on the n-channel-convention quantity, which is
        ## also what SPICE's own `mos1load.c` limits.
        if nmos > 0:
            bgs, bds, bsb = Branch(g, s), Branch(d, s), Branch(s, b)
            bdb, bgb = Branch(d, b), Branch(g, b)
        else:
            bgs, bds, bsb = Branch(s, g), Branch(s, d), Branch(b, s)
            bdb, bgb = Branch(b, d), Branch(b, g)

        ## -- temperature (EKV 2.6, section "Temperature effects") ------
        ut = _var(_vt(T), 'ut')
        trat = _var(T / tnom, 'trat')                              # noqa
        ltr = _var(sympy.log(trat), 'ltrat')
        egT = _var(1.16 - 7.02e-4 * T ** 2 / (T + 1108.0), 'egT')
        egn = _var(1.16 - 7.02e-4 * tnom ** 2 / (tnom + 1108.0),   # noqa
                   'egtnom')
        vtoT = _var(vto - tcv * (T - tnom), 'vtoT')                # noqa
        ## `safe_pow` with the base floored relative to the ratio's own
        ## scale: `trat` is a temperature ratio and cannot legitimately
        ## be below 1e-3 (0.3 K), and `bex` is negative on every real
        ## card, so an unfloored base is a genuine pole.
        kpT = _var(kp * _safe_pow(trat, bex, lo=1e-3), 'kpT')      # noqa
        phiT = _var(phi * trat - 3.0 * ut * ltr                    # noqa
                    - egn * trat + egT, 'phiT')

        ## -- the three probes ------------------------------------------
        ## `limit_fet` on the gate drive and `limit_vds` on the drain,
        ## which is SPICE's pair.  WHICH terminal each write-back lands on
        ## is decided at run time (the one that drifted further from the
        ## last accepted point -- roadmap 12.1), so neither declaration
        ## order nor a fixed node assignment enters into it.
        ##
        ## `von`, the BODY-BIASED turn-on, and not `vto`: SPICE passes
        ## `fetlim` a `von` recomputed from the previous iterate's bulk
        ## bias, and a limiter parameter that reads the solution is now
        ## evaluated at exactly that point -- the last accepted iterate.
        ## Until 2026-08-25 this was refused and the model passed its
        ## zero-bias `vto`, which at 2 V of body bias put every clamp
        ## 565 mV below where the device actually turns on.  `maxc`
        ## inside the root because the parameter is evaluated NUMERICALLY
        ## at a point Newton may have pushed into forward body bias.
        ## THE BULK PROBE.  The current reads `vgs`, `vds` AND `vsb`;
        ## SPICE's limiter set covers the first two, and the declared-
        ## probe route refuses a device whose current reaches the
        ## solution through an unlimited branch -- it said "var(vsb)
        ## reads b, g, which no $limit probe limits".  So the third
        ## branch needs a probe for the model to enter vector PCNR at
        ## all, and with one the device owns all three of its unknowns
        ## and the tail-node clash cannot form on it (roadmap sec. 15).
        ##
        ## It was `limit_identity` -- a probe with NO law -- from Stage 2
        ## until 2026-08-26.  That admitted the model and nothing more,
        ## and sec. 27's wild-start fix is precisely to pass the SEED
        ## through each device's own law: an identity law clamps nothing,
        ## so this was the one circuit left falling back to the ordinary
        ## solver from a wild start.
        ##
        ## `limit_delta` (sec. 34) gives it a real law without inventing
        ## physics the model does not have.  Measured on the diff pair,
        ## seven `vin` values x both instance orders x a uniform 20 V
        ## start: identity fails **14 of 14**, this converges **14 of
        ## 14**, order-independently, with every zero-start result and
        ## every converged tail unchanged.
        ##
        ## ⚠ 10 V IS CHOSEN FROM TWO GRIDS, NOT A POINT -- the bound is
        ## not monotone, and the two circuits that care disagree about
        ## which end of the range is good.
        ##
        ## Diff pair, wild start (14 cases): <= 2 V partial, and
        ## ORDER-DEPENDENT at |vin| = 1, which is the one property PCNR
        ## exists to provide; **3 to 20 V all fourteen**; 50 V fails
        ## everything, being too loose to limit at all.
        ##
        ## 40 V cascode from the origin, Jacobian evaluations (the
        ## `test_fet_limiting_cuts_the_iteration_count_on_a_hard_solve`
        ## benchmark): identity 13, delta 1 -> 20, 2 -> 16, 3 -> 13,
        ## **5 -> 14**, 10 -> 13, 20 -> 12.
        ##
        ## So 5 V -- picked first, from the diff-pair grid alone -- is
        ## the ONE value inside the plateau that costs this solve an
        ## iteration. 10 V is inside both: all fourteen on the pair, and
        ## it ties the identity probe on the cascode.  1 V -- the value
        ## `compact.PspMosLongChannel.limit` uses by hand, and the
        ## obvious guess -- is the worst of the lot at 20.
        vsb = _var(limit_delta(bsb.V, 10.0), 'vsb')
        von = _var(vtoT + gamma * (sympy.sqrt(_maxc(phiT + vsb, 0.0))  # noqa
                                   - sympy.sqrt(phiT)), 'von')
        vgs = _var(limit_fet(bgs.V, von), 'vgs')
        vds = _var(limit_vds(bds.V), 'vds')
        ## Referred to the BULK, which is what makes the model
        ## symmetric in source and drain.
        vgb = _var(vgs + vsb, 'vgb')
        vdb = _var(vds + vsb, 'vdb')

        ## -- pinch-off voltage and slope factor ------------------------
        vgp = _var(vgb - vtoT + phiT + gamma * sympy.sqrt(phiT),   # noqa
                   'vgp')
        ## `maxc(vgp, 0)` inside the root, because the arm below is
        ## selected for `vgp <= 0` and the root's argument goes negative
        ## there.  The `1e-12` is not a smoothing: it exists so that the
        ## IDEAL body-effect-free device (`gamma = 0`, a legitimate card
        ## and the one the construction tests use) does not evaluate
        ## `0 * d(sqrt(0))/dvgp` = `0 * inf` = NaN.  With any real
        ## `gamma` the `gamma^2/4` term dominates it by twenty decades.
        vgr = _var(sympy.sqrt(_maxc(vgp, 0.0) + 0.25 * gamma ** 2  # noqa
                              + 1e-12), 'vgr')
        ## EKV's own `if (VGprime > 0) ... else VP = -PHI`.  The join is
        ## C1: the value is `-PHI` and the slope is `1 - gamma/gamma = 0`
        ## on both sides of `vgp = 0`.
        vp = _var(sympy.Piecewise(
            (vgp - phiT - gamma * (vgr - 0.5 * gamma), vgp > 0.0),  # noqa
            (-phiT, True)), 'vp')
        ## `vp + phi >= 0` STRUCTURALLY (the `vgp > 0` arm is the
        ## increasing function `vgp - gamma*sqrt(vgp + g^2/4) + g^2/2`,
        ## which is zero at `vgp = 0`), and `4*ut` is 0.1 V on top of
        ## that.  So this square root needs no regulariser, and it does
        ## not get one: every `safe_` helper spends exponent range, and
        ## spending it on a guard that cannot fire is a pure loss
        ## (hdl.md 3.2c).
        nq = _var(1.0 + gamma / (2.0 * sympy.sqrt(vp + phiT        # noqa
                                                  + 4.0 * ut)), 'nq')

        ## -- the interpolation function --------------------------------
        beta = _var(kpT * w / l / (1.0 + theta * _maxc(vp, 0.0)),  # noqa
                    'beta')
        ispec = _var(2.0 * nq * beta * ut * ut, 'ispec')
        xf = _var((vp - vsb) / (2.0 * ut), 'xf')
        xr = _var((vp - vdb) / (2.0 * ut), 'xr')
        sf = _var(_softplus(xf), 'sf')
        sr = _var(_softplus(xr), 'sr')
        iff = _var(sf * sf, 'iff')
        irr = _var(sr * sr, 'irr')
        ids = _var(ispec * (iff - irr), 'ids')

        ## -- the four terminal charges ---------------------------------
        ## The charge-based result: with `i = q^2 + q` the normalised
        ## charge density at a channel end is `q = sqrt(1/4 + i) - 1/2`,
        ## and the Ward-Dutton partition of the total inversion charge
        ## integrates in closed form to the two cubics below.  Their
        ## normalisation is `W*L*cox*n*UT`, which is fixed by the
        ## symmetric case: at `Vds = 0` both reduce to `q_s`, and
        ## `Q_I = -W*L*cox*n*UT*(qsn + qdn)` is then the uniform-channel
        ## charge `-n*cox*(VP - VS)*W*L`.
        ##
        ## `0.25 + iff >= 0.25` because `iff` is a square, so these
        ## roots need no guard either, and `(sxf + sxr)^2 >= 1` for the
        ## same reason, so neither does the division.
        sxf = _var(sympy.sqrt(0.25 + iff), 'sxf')
        sxr = _var(sympy.sqrt(0.25 + irr), 'sxr')
        sden = _var((sxf + sxr) ** 2, 'sden')
        qdn = _var(4.0 / 15.0 * (3.0 * sxr ** 3 + 6.0 * sxr ** 2 * sxf
                                 + 4.0 * sxr * sxf ** 2 + 2.0 * sxf ** 3)
                   / sden - 0.5, 'qdn')
        qsn = _var(4.0 / 15.0 * (3.0 * sxf ** 3 + 6.0 * sxf ** 2 * sxr
                                 + 4.0 * sxf * sxr ** 2 + 2.0 * sxr ** 3)
                   / sden - 0.5, 'qsn')
        cwl = _var(w * l * cox, 'cwl')                             # noqa
        qd = _var(-cwl * nq * ut * qdn, 'qd')
        qs = _var(-cwl * nq * ut * qsn, 'qs')
        ## The bulk charge, linearised about pinch-off exactly as the
        ## current is: `dQb'/dV = -cox*gamma/(2*sqrt(psi)) = -cox*(n-1)`
        ## and `V - VP = Qinv'/(n*cox)`, so
        ## `Qb' = -cox*gamma*sqrt(phi + VP) - (n-1)/n*Qinv'`.  The
        ## depletion root IS reachable at zero -- `vp = -phi` is deep
        ## accumulation -- so this one is smoothed, with a width
        ## relative to `phi` rather than absolute.
        qb = _var(-cwl * gamma * sympy.sqrt(                       # noqa
            _hypsmooth(phiT + vp, 1e-6 * phi))                     # noqa
            - (nq - 1.0) / nq * (qd + qs), 'qb')
        ## Charge neutrality is STRUCTURAL: the gate takes what is left,
        ## so the four sum to zero to the last bit whatever the physics
        ## above does.
        qg = _var(-(qd + qs + qb), 'qg')

        ## -- noise -----------------------------------------------------
        ## Channel thermal noise from the Klaassen-Prins integral,
        ## `S = 4*k*T*mu*|Q_I|/L^2`, which in these variables is
        ## `4*k*T*beta*n*UT*(qsn + qdn)`.  Two limits check it: at
        ## `Vds = 0` it is exactly `4*k*T*gds` (Nyquist), and in
        ## saturation it is `(2/3)*4*k*T*n*gm`.  Both are asserted in
        ## the tests rather than the formula being taken on trust.
        gn = _var(beta * nq * ut * _maxc(qsn + qdn, 0.0), 'gn')
        noise = (
            Contribution(bds.I, _white_noise(4.0 * _KB * T * gn)),
            ## EKV normalises flicker noise by `cox*W*L`; with `cox`
            ## defaulting to zero that would be a division by zero, so
            ## the normalisation is folded into `kf` here as it is in
            ## the diode above.  `safe_abs` because a PSD may not be
            ## negative and `ids` is negative for a reversed device.
            Contribution(bds.I, _flicker_noise(
                kf * _safe_abs(ids) ** af, ef)),                   # noqa
        )
        return (Contribution(bds.I, ids),
                Contribution(bdb.I, ddt(qd)),
                Contribution(bgb.I, ddt(qg)),
                Contribution(bsb.I, ddt(qs))) + noise
    return analog


class EkvNmosHdl(Behavioural):
    """EKV 2.6 long-channel n-channel MOSFET.  Terminals ``(d, g, s, b)``.

    One expression covers weak, moderate and strong inversion; see
    `_ekv_analog` for the equations, the references and the list of
    what is deliberately absent.

    This model exists in the library for two reasons beyond being a
    MOSFET.  It is the first production user of `hdl.limit_fet` and
    `hdl.limit_vds`, and of `hdl.softplus`.  And it is an INDEPENDENT
    FORMULATION of the same physics as `compact.PspMosLongChannel`:
    PSP solves the surface potential and EKV linearises the depletion
    charge about pinch-off, so where the two agree the agreement is
    evidence about both, and where they disagree the size of the
    disagreement is the price of the linearisation.  The cross-check is
    in ``test_elements_hdl_library3.py`` and its result is reported
    there rather than asserted here.
    """
    instparams = _ekv_params()

    analog = staticmethod(_ekv_analog(TEMP, +1))


class EkvPmosHdl(Behavioural):
    """EKV 2.6 long-channel p-channel MOSFET.  Terminals ``(d, g, s, b)``.

    `EkvNmosHdl` with every branch declared with its terminals
    exchanged; the card carries positive magnitudes for both devices
    (``vto`` is the magnitude of the threshold, ``gamma`` and ``phi``
    are positive) and the device conducts for
    ``V(s) > V(g) + vto`` and ``V(s) > V(d)``.

    Redefining ``analog`` here is REQUIRED and not decorative -- see
    `GummelPoonPnpHdl`.
    """
    instparams = _ekv_params()

    analog = staticmethod(_ekv_analog(TEMP, -1))


## ======================================================================
## SPICE MOS level 1 (Shichman-Hodges) -- roadmap item 10
## ======================================================================

#: Permittivities and the intrinsic carrier concentration, SI, as SPICE
#: uses them to derive `gamma` and `phi` from `tox` and `nsub`.  `nsub`
#: is a card value in cm^-3 -- every SPICE card in existence writes it
#: that way -- so it is converted here rather than at the call site.
_EPSOX = 3.9 * 8.854187817e-12
_EPSSI = 11.7 * 8.854187817e-12
_NI_CM3 = 1.45e10


def _mos1_params():
    """The SPICE level-1 MOSFET card, in SPICE's own names and defaults.

    **Three of SPICE's card names cannot be Python identifiers**, and the
    DSL's bare-name convention makes a parameter name exactly that -- the
    compiler injects the symbols into `analog()`'s globals, so a model
    reads a parameter by writing its name.  So:

    ==============  ==========  ==========================================
    SPICE           here        why
    ==============  ==========  ==========================================
    ``IS``          ``IS``      ``is`` is a keyword (`DiodeSpiceHdl`'s
                                spelling, kept)
    ``LAMBDA``      ``lambd``   ``lambda`` is a keyword
    ``AS``          ``asrc``    ``as`` is a keyword
    ==============  ==========  ==========================================

    ``globals()['lambda']`` would in fact resolve inside ``analog()``, and
    ``cls(d, g, s, b, **{'lambda': 0.02})`` would set it, so the name is
    reachable -- but neither spelling is one anybody should have to write,
    and a parameter that cannot be typed is worse than a renamed one.
    Recorded in the batch's friction log.  Since 2026-08-26 SPICE's own
    spellings are accepted on the instance as ALIASES
    (``aliasparams = {'lambda': 'lambd', 'as': 'asrc'}`` on both
    classes, `Behavioural.aliasparams`), so a card can be passed as
    written; the canonical names stay ``lambd``/``asrc`` because this
    model reads its parameters as bare names.  A ``params_as`` model
    could declare ``lambda`` outright and read it as ``p['lambda']``
    (`hdl.ParamNamespace`).

    ``gamma`` and ``phi`` follow SPICE's rule: **given on the card they are
    used; absent they are derived** from ``nsub`` and ``tox``.  That is
    what `hdl.param_given` is for, and this is its first production user
    alongside ``$limit`` -- a combination that was a run-time
    ``TypeError`` until 2026-08-25 (`test_limit_with_param_given.py`).

    ``ad``/``asrc``/``pd``/``ps``/``nrd``/``nrs`` are SPICE's per-instance
    geometry, and without them ``cj``, ``cjsw``, ``js`` and ``rsh`` are
    dead knobs: each of those four is a *density*, and a density with no
    area is nothing.
    """
    return [
        Parameter(name='vto', desc='Zero-bias threshold voltage', unit='V',
                  default=0.0),
        Parameter(name='kp', desc='Transconductance parameter mu*Cox',
                  unit='A/V^2', default=2e-5),
        Parameter(name='gamma', desc='Body effect factor (unset: from '
                  'nsub and tox)', unit='V^0.5', default=0.0),
        Parameter(name='phi', desc='Surface inversion potential (unset: '
                  'from nsub)', unit='V', default=0.6),
        Parameter(name='lambd', desc="Channel-length modulation (SPICE's "
                  'LAMBDA)', unit='1/V', default=0.0),
        Parameter(name='tox', desc='Oxide thickness', unit='m',
                  default=1e-7),
        Parameter(name='nsub', desc='Substrate doping', unit='cm^-3',
                  default=0.0),
        Parameter(name='w', desc='Channel width', unit='m', default=1e-4),
        Parameter(name='l', desc='Drawn channel length', unit='m',
                  default=1e-4),
        Parameter(name='ld', desc='Lateral diffusion', unit='m',
                  default=0.0),
        Parameter(name='cgso', desc='Gate-source overlap capacitance per '
                  'metre of width', unit='F/m', default=0.0),
        Parameter(name='cgdo', desc='Gate-drain overlap capacitance per '
                  'metre of width', unit='F/m', default=0.0),
        Parameter(name='cgbo', desc='Gate-bulk overlap capacitance per '
                  'metre of length', unit='F/m', default=0.0),
        Parameter(name='cbd', desc='Zero-bias bulk-drain capacitance '
                  '(0: use cj*ad)', unit='F', default=0.0),
        Parameter(name='cbs', desc='Zero-bias bulk-source capacitance '
                  '(0: use cj*asrc)', unit='F', default=0.0),
        Parameter(name='IS', desc='Bulk junction saturation current '
                  "(SPICE's IS; used when js*area is zero)", unit='A',
                  default=1e-14),
        Parameter(name='pb', desc='Bulk junction potential at tnom',
                  unit='V', default=0.8),
        Parameter(name='cj', desc='Bulk junction bottom capacitance per '
                  'area', unit='F/m^2', default=0.0),
        Parameter(name='cjsw', desc='Bulk junction sidewall capacitance '
                  'per metre', unit='F/m', default=0.0),
        Parameter(name='mj', desc='Bottom grading coefficient', unit='',
                  default=0.5),
        Parameter(name='mjsw', desc='Sidewall grading coefficient',
                  unit='', default=0.33),
        Parameter(name='fc', desc='Forward-bias depletion coefficient',
                  unit='', default=0.5),
        Parameter(name='js', desc='Bulk junction saturation current per '
                  'area', unit='A/m^2', default=0.0),
        Parameter(name='rd', desc='Drain ohmic resistance (0: use '
                  'rsh*nrd)', unit='ohm', default=0.0),
        Parameter(name='rs', desc='Source ohmic resistance (0: use '
                  'rsh*nrs)', unit='ohm', default=0.0),
        Parameter(name='rsh', desc='Drain/source sheet resistance',
                  unit='ohm/sq', default=0.0),
        Parameter(name='ad', desc='Drain diffusion area', unit='m^2',
                  default=0.0),
        Parameter(name='asrc', desc="Source diffusion area (SPICE's AS)",
                  unit='m^2', default=0.0),
        Parameter(name='pd', desc='Drain diffusion perimeter', unit='m',
                  default=0.0),
        Parameter(name='ps', desc='Source diffusion perimeter', unit='m',
                  default=0.0),
        Parameter(name='nrd', desc='Drain squares', unit='', default=1.0),
        Parameter(name='nrs', desc='Source squares', unit='', default=1.0),
        Parameter(name='kf', desc='Flicker-noise coefficient', unit='',
                  default=0.0),
        Parameter(name='af', desc='Flicker-noise exponent', unit='',
                  default=1.0),
        Parameter(name='tnom', desc='Parameter measurement temperature',
                  unit='K', default=300.15),
    ]


def _mos1_analog(T, nmos, limiting='group'):
    """Build the ``analog()`` body of a SPICE level-1 MOSFET.

    Equations: Massobrio & Antognetti, *Semiconductor Device Modeling
    with SPICE*, 2nd ed., ch. 3 -- Shichman and Hodges, *IEEE J.
    Solid-State Circuits* **3**, 285 (1968), as SPICE2 packaged it --
    together with ngspice's `mos1temp.c` for the temperature path.

    **This model exists in the library to be the first production user of
    `hdl.limit_together`.**  SPICE's own MOSFET limiting is four probes --
    ``fetlim(vgs)``, ``limvds(vds)``, ``pnjlim(vbs)``, ``pnjlim(vbd)`` --
    and a four-probe device **cannot be declared per-probe at all**: each
    per-probe write-back claims a terminal of its own, the four probes
    span four terminals with a cycle in them
    (``(b,s)-(b,d)-(d,s)``), and `generate_code` refuses.  Grouped, the
    write-back is a spanning forest over the terminals and the model
    compiles.  `EkvNmosHdl` declares two probes and is the reason the
    per-probe form exists; this one declares SPICE's four.

    **The source/drain swap, written without a branch.**  SPICE decides
    which diffusion is the source by comparing ``vds`` against zero and
    swaps the terminals, which a compiled expression cannot do (the
    condition is a solution variable).  The usual substitute is
    ``sign(vds) * f(|vds|)``, and it is WRONG HERE: `hdl.sign` has
    ``fdiff = 0`` and ``sign(0) = 0``, so the analytic ``dId/dvds`` at
    ``vds = 0`` -- a switch held fully on, the single most common bias in
    a digital netlist -- comes out as exactly zero where it should be
    ``beta*(vgs - vth)``.  Level 1 does not need it.  Written as

        Id = (beta/2) * [ max(vgs - vth, 0)^2 - max(vgd - vth, 0)^2 ]

    with **one** threshold ``vth``, taken at the more forward-biased of
    the two bulk junctions (which is SPICE's "the source is the lower
    terminal"), the expression is already antisymmetric under exchanging
    source and drain, because exchanging them exchanges the two squares.
    The bias-dependent square root appears once, so no ``|vds|`` and no
    ``sign`` is needed, and the Jacobian at ``vds = 0`` is right.
    ``test_level1_gds_at_vds_zero_is_the_on_conductance`` asserts exactly
    that, against the ``sign``-based form which it also builds.

    The one place ``|vds|`` survives is the channel-length modulation
    factor ``1 + lambda*|vds|``, and there it is harmless: it multiplies
    ``Id``, which is zero wherever ``|vds|`` is non-differentiable, so the
    product rule kills the offending term.

    **Implemented:** the Shichman-Hodges current with body effect, channel
    -length modulation and lateral diffusion; ``gamma``/``phi`` derived
    from ``nsub``/``tox`` when the card does not give them; both bulk
    junctions' current and depletion charge (bottom and sidewall, with
    the ``fc`` linearisation) from either the absolute ``cbd``/``cbs`` or
    the ``cj``/``cjsw`` densities; drain and source ohmic resistance from
    either ``rd``/``rs`` or ``rsh``*squares, on collapsible internal
    nodes; the three overlap capacitances; SPICE's temperature path;
    channel thermal noise from the Klaassen-Prins integral, flicker noise
    and the parasitics' thermal noise.

    **Not implemented, named rather than left to be discovered:** the
    INTRINSIC gate charge.  SPICE level 1 models it with the Meyer
    capacitances, and Meyer is a capacitance model, not a charge model --
    it is famously non-charge-conserving, and this DSL contributes
    ``ddt(q)``, so transcribing Meyer here would mean integrating three
    capacitances that are not the gradient of any charge.  The overlap
    capacitances ARE charges and are included; a card that needs the
    intrinsic gate capacitance wants `EkvNmosHdl`, whose charge model is
    a genuine Ward-Dutton partition.  Also absent: ``uo`` (so ``kp`` is
    never derived from mobility), the substrate current, and ``theta``
    (not a level-1 parameter).
    """
    def analog(d, g, s, b):
        ## Two internal nodes, one per ohmic resistance, each collapsed
        ## away when its resistance is zero -- SPICE's default for both,
        ## so a bare instance has no internal nodes at all.
        di, si = Node('di'), Node('si')
        brd, brs = Branch(d, di), Branch(s, si)

        ## POLARITY LIVES IN THE BRANCH TERMINAL ORDER, exactly as
        ## `_gummel_poon` and `_ekv_analog` do it: a p-channel's branches
        ## are declared reversed, so every probe below is the
        ## n-channel-convention quantity and the whole body is written
        ## once.  It is not only tidiness -- `limit_fet`, `limit_vds` and
        ## `limit_pnj` bound a branch potential IN THE SOLUTION VECTOR, so
        ## a p-channel whose gate drive is `-V(g,s)` would have every one
        ## of the four laws compress the wrong direction.
        if nmos > 0:
            bgs, bds = Branch(g, si), Branch(di, si)
            bbs, bbd = Branch(b, si), Branch(b, di)
            bgd, bgb = Branch(g, di), Branch(g, b)
        else:
            bgs, bds = Branch(si, g), Branch(si, di)
            bbs, bbd = Branch(si, b), Branch(di, b)
            bgd, bgb = Branch(di, g), Branch(b, g)

        ## -- geometry ---------------------------------------------------
        ## `l - 2*ld` is floored: `ld` is a card value and a card with
        ## `2*ld >= l` is a card for a device that does not exist, but the
        ## model divides by the difference.  1 nm is four decades below
        ## any level-1 geometry.
        leff = _var(_maxc(l - 2.0 * ld, 1e-9), 'leff')             # noqa
        ## `tox` is a divisor twice (the derived `gamma`, the flicker
        ## normalisation), so it is floored at a picometre -- below one
        ## silicon lattice spacing, i.e. unreachable by a real card and
        ## finite for `tox = 0`.
        cox = _var(_EPSOX / _maxc(tox, 1e-12), 'cox')              # noqa

        ## -- gamma and phi: given, or derived from nsub and tox ---------
        ## SPICE's rule, and the reason `param_given` is here.  BOTH arms
        ## are evaluated whatever the card says, so the derived forms must
        ## be finite at `nsub = 0` (the default): the root of zero is
        ## zero, and the logarithm is floored at 1 so it returns zero
        ## rather than -inf.
        vtn = _var(_KB * tnom / _QE, 'vtnom')                      # noqa
        nsm3 = _var(_maxc(nsub, 0.0) * 1e6, 'nsm3')                # noqa
        gd = _var(sympy.sqrt(2.0 * _QE * _EPSSI * nsm3) / cox, 'gamd')
        phid = _var(2.0 * vtn * sympy.log(_maxc(nsm3 / (_NI_CM3 * 1e6),
                                                1.0)), 'phid')
        gg = _var(param_given('gamma'), 'ggam')
        gp = _var(param_given('phi'), 'gphi')
        gam = _var(gg * gamma + (1.0 - gg) * gd, 'gam')            # noqa
        ## SPICE clamps PHI to 0.1 V and so does this: the model divides
        ## by `phi` in the forward-bulk continuation and takes its square
        ## root, and a card that says `phi = 0` (or a derived zero, which
        ## `nsub = 0` gives) is not a card for a MOSFET.
        phz = _var(_maxc(gp * phi + (1.0 - gp) * phid, 0.1), 'phz')  # noqa

        ## -- temperature (ngspice `mos1temp.c`) ------------------------
        vtT = _var(_vt(T), 'vtT')
        trat = _var(T / tnom, 'trat')                              # noqa
        ltr = _var(sympy.log(trat), 'ltrat')
        egT = _var(1.16 - 7.02e-4 * T ** 2 / (T + 1108.0), 'egT')
        egn = _var(1.16 - 7.02e-4 * tnom ** 2 / (tnom + 1108.0),   # noqa
                   'egtnom')
        ## The built-in potentials move exactly as the diode's `vjT` does.
        phiT = _var(_maxc(phz * trat - 3.0 * vtT * ltr
                          - egn * trat + egT, 0.1), 'phiT')
        pbT = _var(pb * trat - 3.0 * vtT * ltr                     # noqa
                   - egn * trat + egT, 'pbT')
        ## `vto(T) = vbi(T) + gamma*sqrt(phi(T))` with
        ## `vbi(T) = vbi + (phi(T) - phi)/2`, which is `mos1temp.c`.
        vtoT = _var(vto + 0.5 * (phiT - phz)                       # noqa
                    + gam * (sympy.sqrt(phiT) - sympy.sqrt(phz)), 'vtoT')
        ## Mobility falls as T^-1.5, so `kp` does.  `safe_pow` with the
        ## base floored relative to the ratio's own scale -- a temperature
        ## ratio cannot legitimately be below 1e-3 (0.3 K) and the
        ## exponent is negative, so an unfloored base is a genuine pole.
        kpT = _var(kp * _safe_pow(trat, -1.5, lo=1e-3), 'kpT')     # noqa
        ## ngspice's own junction-current temperature law for the MOS
        ## (`mos1temp.c`): `IS(T) = IS*exp(-Eg(T)/Vt(T) + Eg(Tnom)/Vt(Tnom))`,
        ## which is exactly 1 at `T = tnom` and needs no `xti` -- level 1
        ## has no such parameter.
        factlog = _var(-egT / vtT + egn / vtn, 'factlog')
        ## The junction capacitances' temperature factors, the diode's,
        ## one per grading coefficient because SPICE scales the bottom
        ## and the sidewall with their own.
        tshift = _var(4e-4 * (T - tnom) - (pbT / pb - 1.0), 'tshift')  # noqa
        cfacb = _var(1.0 + mj * tshift, 'cfacb')                   # noqa
        cfacw = _var(1.0 + mjsw * tshift, 'cfacw')                 # noqa

        ## -- the junction saturation currents --------------------------
        ## SPICE's rule: the per-area `js` if the card gives an area for
        ## it, the absolute `IS` otherwise.  A parameter-only `Piecewise`,
        ## not `param_given`, because the condition SPICE actually applies
        ## is "is the PRODUCT non-zero", which a givenness flag cannot
        ## express.  Floored so that `limit_pnj` can take its logarithm
        ## on a card that switches both off.
        isbd = _var(_maxc(sympy.Piecewise((js * ad, js * ad > 0.0),  # noqa
                                          (IS, True))              # noqa
                          * _expl(factlog), 1e-30), 'isbd')
        isbs = _var(_maxc(sympy.Piecewise((js * asrc, js * asrc > 0.0),
                                          (IS, True))              # noqa
                          * _expl(factlog), 1e-30), 'isbs')

        ## -- the four probes, LIMITED AS ONE ---------------------------
        ## SPICE's own set.  Four probes over four terminals, and
        ## `(b,s)`, `(b,d)` and `(d,s)` close a triangle, so the per-probe
        ## write-back has no terminal left for the fourth and
        ## `generate_code` refuses the model outright.  `limit_together`
        ## is what makes the declaration expressible; see `hdl.
        ## limit_together` and roadmap 10.3(b).
        ##
        ## `von`, the BODY-BIASED turn-on SPICE passes to `fetlim`, built
        ## from the RAW bulk-source potential with the same continuation
        ## `vth` below uses.  A limiter parameter that reads the solution
        ## is evaluated at the last accepted iterate (roadmap 12.6(c),
        ## closed 2026-08-25), which is exactly SPICE's `von` semantics;
        ## until then this passed `vtoT` and was measurably loose.
        ##
        ## `limiting` selects the declaration, and the three values exist
        ## so that the limiter can be measured against its ABSENCE on
        ## THIS model rather than on a toy (`newton-limiting`: "measure
        ## the intervention against its absence, on a real solve, at the
        ## scale it is used").  `'group'` is what the shipped classes
        ## use; `'probe'` is the per-probe declaration, which does not
        ## compile and is kept so that a test can say so; `'none'` is the
        ## unlimited device.
        vbsl = _var(bbs.V, 'vbsl')
        sargl = _var(_select(
            (sympy.sqrt(phiT - vbsl), vbsl <= 0.0),
            (sympy.sqrt(phiT) / (1.0 + vbsl / (2.0 * phiT)),
             True)), 'sargl')
        von = _var(vtoT + gam * (sargl - sympy.sqrt(phiT)), 'von')
        if limiting == 'group':
            vgs, vds, vbs, vbd = limit_together(                   # noqa
                limit_fet(bgs.V, von), limit_vds(bds.V),
                limit_pnj(bbs.V, isbs, vtT), limit_pnj(bbd.V, isbd, vtT))
        elif limiting == 'probe':
            vgs, vds = limit_fet(bgs.V, von), limit_vds(bds.V)
            vbs = limit_pnj(bbs.V, isbs, vtT)
            vbd = limit_pnj(bbd.V, isbd, vtT)
        else:
            vgs, vds, vbs, vbd = bgs.V, bds.V, bbs.V, bbd.V
        vgs = _var(vgs, 'vgs')
        vds = _var(vds, 'vds')
        vbs = _var(vbs, 'vbs')
        vbd = _var(vbd, 'vbd')
        vgd = _var(vgs - vds, 'vgd')

        ## -- the threshold ---------------------------------------------
        ## Taken at the MORE FORWARD-BIASED junction, which is SPICE's
        ## "the source is the lower diffusion" written without a branch on
        ## the solution.
        vbx = _var(_maxc(vbs, vbd), 'vbx')
        ## SPICE's own continuation for a forward-biased bulk
        ## (`mos1load.c`): the square root below zero bias, its first-order
        ## expansion above it.  Both arms have to be finite in the region
        ## they are NOT selected in, and `hdl.select` derives that from
        ## the condition rather than the author choosing a floor: inside
        ## the first arm `vbx` becomes `minc(vbx, 0)`, so the radicand is
        ## at least `phiT`; inside the second it becomes `maxc(vbx, 0)`,
        ## so the denominator is at least 1.  Where each arm IS selected
        ## the clamp is the identity -- bit for bit, and in the
        ## derivative too -- so this is the same model.
        ##
        ## What it replaces: two hand-picked constants (`0.25*phiT` and
        ## `0.1`) which were correct, provably never reached, and
        ## unrelated to the bound the condition already states.  The
        ## same four lines appear at four sites (this one, `sargl`
        ## above, and both of `_mos3_channel`/`_mos3_analog`).
        sarg = _var(_select(
            (sympy.sqrt(phiT - vbx), vbx <= 0.0),
            (sympy.sqrt(phiT) / (1.0 + vbx / (2.0 * phiT)),
             True)), 'sarg')
        vth = _var(vtoT + gam * (sarg - sympy.sqrt(phiT)), 'vth')
        vgt = _var(_maxc(vgs - vth, 0.0), 'vgt')
        vgtd = _var(_maxc(vgd - vth, 0.0), 'vgtd')

        ## -- the current ------------------------------------------------
        beta = _var(kpT * w / leff, 'beta')                        # noqa
        ## `|vds|` written as a difference of the two bulk-referred
        ## potentials, so it is one expression rather than a `sign`.  It
        ## multiplies `ids0`, which is zero wherever this is kinked.
        vdsa = _var(_maxc(vbs, vbd) - _minc(vbs, vbd), 'vdsa')
        clm = _var(1.0 + lambd * vdsa, 'clm')                      # noqa
        ids0 = _var(0.5 * beta * (vgt * vgt - vgtd * vgtd), 'ids0')
        ids = _var(ids0 * clm, 'ids')

        ## -- the two bulk junctions -------------------------------------
        ibs = _var(isbs * (_expl(vbs / vtT) - 1.0), 'ibs')
        ibd = _var(isbd * (_expl(vbd / vtT) - 1.0), 'ibd')
        ## Bottom capacitance: the absolute `cbd`/`cbs` if the card gives
        ## it, `cj` times the area otherwise -- SPICE's rule again.
        cbdb = _var(cfacb * sympy.Piecewise((cbd, cbd > 0.0),       # noqa
                                           (cj * ad, True)), 'cbdb')  # noqa
        cbsb = _var(cfacb * sympy.Piecewise((cbs, cbs > 0.0),       # noqa
                                           (cj * asrc, True)), 'cbsb')
        cbdw = _var(cfacw * cjsw * pd, 'cbdw')                      # noqa
        cbsw = _var(cfacw * cjsw * ps, 'cbsw')                      # noqa
        qbd = _var(_pn_depletion_charge(vbd, cbdb, pbT, mj, fc, 'bd')
                   + _pn_depletion_charge(vbd, cbdw, pbT, mjsw, fc,  # noqa
                                          'bdw'), 'qbd')           # noqa
        qbs = _var(_pn_depletion_charge(vbs, cbsb, pbT, mj, fc, 'bs')
                   + _pn_depletion_charge(vbs, cbsw, pbT, mjsw, fc,  # noqa
                                          'bsw'), 'qbs')           # noqa

        ## -- the ohmic resistances --------------------------------------
        ## `rd` if given, `rsh` times the squares otherwise; the branch
        ## collapses only when BOTH are zero, which is SPICE's default.
        rdx = _var(sympy.Piecewise((rd, rd > 0.0),                 # noqa
                                   (rsh * nrd, True)), 'rdx')      # noqa
        rsx = _var(sympy.Piecewise((rs, rs > 0.0),                 # noqa
                                   (rsh * nrs, True)), 'rsx')      # noqa

        ## -- statements -------------------------------------------------
        stmts = (
            Contribution(bds.I, ids),
            Contribution(bbs.I, ibs + ddt(qbs)),
            Contribution(bbd.I, ibd + ddt(qbd)),
            ## The three overlap capacitances.  These are CHARGES, which
            ## is why they are here and the intrinsic Meyer capacitances
            ## are not: `cgso*w*V` is the gradient of `cgso*w*V^2/2` and
            ## Meyer's three are the gradient of nothing.
            Contribution(bgs.I, ddt(cgso * w * bgs.V)),            # noqa
            Contribution(bgd.I, ddt(cgdo * w * bgd.V)),            # noqa
            Contribution(bgb.I, ddt(cgbo * leff * bgb.V)),         # noqa
            Contribution(brd.I, brd.V / rdx),
            Collapse(brd, sympy.And(rd <= 0.0, rsh * nrd <= 0.0)),  # noqa
            Contribution(brs.I, brs.V / rsx),
            Collapse(brs, sympy.And(rs <= 0.0, rsh * nrs <= 0.0)),  # noqa
        )

        ## -- noise -------------------------------------------------------
        ## Channel thermal noise from the Klaassen-Prins integral for a
        ## square-law channel,
        ##
        ##     S = (8/3)*k*T*beta*(vgt^2 + vgt*vgtd + vgtd^2)/(vgt + vgtd)
        ##
        ## which has both textbook limits built in and neither pasted on:
        ## at `Vds = 0` (vgt = vgtd) it is `4*k*T*beta*vgt`, which is
        ## `4*k*T*gds` exactly (Nyquist), and in saturation (vgtd = 0) it
        ## is `(2/3)*4*k*T*beta*vgt`, which is `(2/3)*4*k*T*gm` exactly.
        ## Both are asserted in the tests rather than the formula being
        ## taken on trust.  `safe_div` because a device in cutoff has
        ## `vgt = vgtd = 0` and the ratio is 0/0 there, where the noise is
        ## zero anyway.
        gn = _var(_safe_div(vgt * vgt + vgt * vgtd + vgtd * vgtd,
                            vgt + vgtd), 'gn')
        noise = (
            Contribution(bds.I, _white_noise(8.0 / 3.0 * _KB * T
                                             * beta * gn)),
            ## SPICE's flicker normalisation, `kf*Id^af/(cox*Leff^2*f)`.
            ## `safe_abs` because a PSD may not be negative and `ids` is
            ## negative for a reversed device.
            Contribution(bds.I, _flicker_noise(
                kf * _safe_abs(ids) ** af / (cox * leff ** 2), 1)),  # noqa
            Contribution(brd.I, _white_noise(4.0 * _KB * T / rdx)),
            Contribution(brs.I, _white_noise(4.0 * _KB * T / rsx)),
        )
        return stmts + noise
    return analog


class MosLevel1Hdl(Behavioural):
    """SPICE level-1 (Shichman-Hodges) n-channel MOSFET.
    Terminals ``(d, g, s, b)``.

    Parameters and defaults are SPICE's -- so a bare
    ``MosLevel1Hdl('d', 'g', 's', 'b')`` is a ``vto = 0``,
    ``kp = 2e-5``, 100 um by 100 um device with no body effect, no
    channel-length modulation, no parasitic resistance and no charge
    except the two bulk junctions' current, and every block is switched
    on by giving its parameter.  See `_mos1_analog` for the equations,
    the reference and what is deliberately absent (the intrinsic Meyer
    capacitances, and why).

    **The first production user of `hdl.limit_together`.**  This model
    declares SPICE's own four limiter probes -- ``fetlim(vgs)``,
    ``limvds(vds)``, ``pnjlim(vbs)``, ``pnjlim(vbd)`` -- which is the
    declaration a per-probe limiter *cannot compile*: four probes over
    four terminals with the triangle ``(b,s)-(b,d)-(d,s)`` in them, so
    the fourth probe finds both of its terminals already claimed by an
    earlier one's write-back and `generate_code` refuses the model.
    `EkvNmosHdl` declares two and is why the per-probe form exists; this
    is the four-probe case the grouped form was built for.  ``explain()``
    prints the group.

    **And it is the first model to use ``$param_given`` and ``$limit``
    together**, which was a run-time ``TypeError`` in generated code until
    2026-08-25: ``gamma`` and ``phi`` follow SPICE's "derived from
    ``nsub`` and ``tox`` unless the card gives them" rule, which is what
    givenness is for and what a default VALUE cannot express.
    """
    instparams = _mos1_params()
    #: SPICE's own spellings of the two keyword-named parameters.
    aliasparams = {'lambda': 'lambd', 'as': 'asrc'}

    analog = staticmethod(_mos1_analog(TEMP, +1))


class MosLevel1PmosHdl(Behavioural):
    """SPICE level-1 p-channel MOSFET.  Terminals ``(d, g, s, b)``.

    `MosLevel1Hdl` with every branch declared with its terminals
    exchanged; the card carries positive magnitudes for both devices
    (``vto`` is the magnitude of the threshold, ``gamma`` and ``phi`` are
    positive) and the device conducts for ``V(s) > V(g) + vto`` and
    ``V(s) > V(d)``.  The two bulk junctions are the n-well's, so they
    conduct for ``V(s) > V(b)``, and their probes are declared that way
    round for exactly the reason `GummelPoonPnpHdl` reverses its --
    ``pnjlim`` must see the voltage that has the exponential in it.

    Redefining ``analog`` here is REQUIRED and not decorative -- see
    `GummelPoonPnpHdl`.
    """
    instparams = _mos1_params()
    aliasparams = {'lambda': 'lambd', 'as': 'asrc'}

    analog = staticmethod(_mos1_analog(TEMP, -1))


## ======================================================================
## Boyle-class operational amplifier macromodel -- roadmap item 4
## ======================================================================


def _opamp_params():
    """The opamp card, in DATASHEET quantities.

    Boyle's own model is parameterised by its internal elements -- a tail
    current, a compensation capacitance, two transistors' betas -- and
    the paper's contribution is the set of closed-form relations that get
    those from the datasheet.  Those relations are what this model keeps
    (see `OpAmpHdl`); the card is what a designer has in front of them.

    ``cc`` is the one internal quantity kept as a parameter, because
    Boyle's relations fix only the two RATIOS ``gbw = gm/(2*pi*cc)`` and
    ``sr = i_tail/cc`` -- so ``cc`` chooses the scale of the gain node and
    nothing observable at the pins depends on it.  It is here so that
    ``sr`` and ``gbw`` can be given independently, which is what a
    datasheet does.
    """
    return [
        Parameter(name='aol', desc='DC open-loop differential gain',
                  unit='V/V', default=2e5),
        Parameter(name='gbw', desc='Gain-bandwidth product', unit='Hz',
                  default=1e6),
        Parameter(name='sr', desc='Slew rate', unit='V/s', default=5e5),
        Parameter(name='cc', desc='Compensation capacitance (sets the '
                  'gain node scale only)', unit='F', default=30e-12),
        Parameter(name='rin', desc='Differential input resistance',
                  unit='ohm', default=2e6),
        Parameter(name='ricm', desc='Common-mode input resistance',
                  unit='ohm', default=1e9),
        Parameter(name='ib', desc='Input bias current', unit='A',
                  default=80e-9),
        Parameter(name='ios', desc='Input offset current', unit='A',
                  default=20e-9),
        Parameter(name='vos', desc='Input offset voltage', unit='V',
                  default=1e-3),
        Parameter(name='cmrr', desc='Common-mode rejection ratio',
                  unit='V/V', default=1e5),
        Parameter(name='psrr', desc='Power-supply rejection ratio',
                  unit='V/V', default=1e5),
        Parameter(name='vsupnom', desc='Total supply voltage at which '
                  'vos is specified', unit='V', default=30.0),
        Parameter(name='rout', desc='Open-loop output resistance',
                  unit='ohm', default=75.0),
        Parameter(name='isc', desc='Short-circuit output current',
                  unit='A', default=25e-3),
        Parameter(name='vdrop', desc='Output saturation voltage below '
                  'each rail', unit='V', default=1.5),
    ]


class OpAmpHdl(Behavioural):
    """A Boyle-class operational amplifier macromodel.
    Terminals ``(inp, inn, out, vcc, vee)``.

    Reference: G. R. Boyle, B. M. Cohn, D. O. Pederson and J. E. Solomon,
    "Macromodeling of integrated circuit operational amplifiers", *IEEE
    J. Solid-State Circuits* **SC-9**, 353-364 (1974).  What survives
    from that paper here is not its topology -- Boyle's input stage is a
    real differential pair, because in 1974 a simulator had transistors
    and nothing else -- but its **relations**, which are what make a
    macromodel a model rather than a curve fit, and which are what the
    tests assert:

    ===================  ================================================
    Boyle's relation     here
    ===================  ================================================
    ``GBW = gm/(2*pi*C)``  ``gm1 = 2*pi*gbw*cc``
    ``SR = I/C``           the gain node is charged through
                           ``isl*tanh(gm1*verr/isl)`` with
                           ``isl = sr*cc``, so no bias can push more
                           than ``sr*cc`` into ``cc``
    ``f_dominant =
      GBW/A_ol``           ``r2 = aol/gm1``, so the gain node's pole is
                           ``1/(2*pi*r2*cc) = gbw/aol`` identically
    ``A_ol = gm*R``        DC gain is ``gm1*r2`` by construction
    ===================  ================================================

    Every one of those is an IDENTITY of the parameterisation, not a fit,
    which is why the tests can assert them to a part in 1e4 off a real AC
    and transient solve rather than to a band.

    **What is here.**  Input stage: differential and common-mode input
    resistance, input bias current, input offset current, input offset
    voltage.  Error amplifier: the differential input plus the
    common-mode error ``vcm/cmrr`` plus the supply error
    ``(vsup - vsupnom)/psrr``, all referred to the input, which is the
    definition of CMRR and PSRR and makes them measurable as ratios of
    two gains.  Gain stage: one dominant pole and slew-rate limiting.
    Output stage: rail limiting a ``vdrop`` short of each supply, output
    resistance and short-circuit current limiting -- **the last two from
    one expression**, ``isc*tanh(dv/(rout*isc))``, whose small-signal
    slope is ``1/rout`` and whose asymptote is ``isc``.

    **What is not, named rather than left to be discovered.**  There is
    **one pole**.  Boyle's model has a second one from the input stage,
    and a designer who wants a phase margin needs it; it is absent here
    because a second pole is a second internal node and the roadmap's
    item is the input/gain/output chain.  Also absent: input-stage
    common-mode range limiting, the asymmetric slew rate real amplifiers
    have (``SR+`` != ``SR-``), input noise, and the supply current, of
    which only the output current and the input bias currents are
    modelled.

    **Two structural choices worth knowing before connecting one.**

    1. *There is no ground pin*, and the output is referenced to the
       **midpoint of the supplies**.  An HDL element can only reach the
       nodes in its own signature, so "referred to the global ground" is
       not expressible from inside ``analog()`` (`SelfHeating` records
       the same limitation for its ambient pin).  For the usual
       ``+-15 V`` supplies the midpoint IS ground and nothing is lost;
       for a single supply it centres the output, which is what a
       single-supply macromodel should do.
    2. *All output current returns to* ``vee``.  A real amplifier's
       output current returns through whichever rail is sourcing it, so
       the supply currents here are right in magnitude and wrong in
       which rail they appear on.  It matters for a supply-current
       measurement and for nothing else.

    ``macromodels.OpAmp`` is the hand-written predecessor; it is a
    four-terminal ``SubCircuit`` with a differential output and no supply
    pins, so this is **not** a drop-in replacement for it.  See
    ``test_elements_hdl_library4.py`` section 2 for what it does and does
    not carry over, measured.
    """
    instparams = _opamp_params()

    @staticmethod
    def analog(inp, inn, out, vcc, vee):
        ## The supply, and the midpoint everything is referred to.
        bsup = Branch(vcc, vee)
        vsup = _var(bsup.V, 'vsup')
        vmid = _var(0.5 * vsup, 'vmid')

        ## -- input stage ------------------------------------------------
        bd = Branch(inp, inn)
        bp, bn = Branch(inp, vee), Branch(inn, vee)
        ## The common-mode input voltage, referred to the supply midpoint.
        vcm = _var(0.5 * (bp.V + bn.V) - vmid, 'vcm')
        ## The three rejection ratios are DIVISORS, so each is floored at
        ## the value below which it stops being a rejection ratio at all
        ## -- 1, i.e. 0 dB.  This is the `differentiable-numerics` rule
        ## in its parameter form: `cmrr = 0` is a card somebody will
        ## write meaning "perfect", and it is a division by zero.  A
        ## floor of 1 turns it into "no rejection", which is at least a
        ## number; the value that means "perfect" is a large one.
        cmr = _var(_maxc(cmrr, 1.0), 'cmrrc')                      # noqa
        psr = _var(_maxc(psrr, 1.0), 'psrrc')                      # noqa
        aolc = _var(_maxc(aol, 1.0), 'aolc')                       # noqa
        ## The error voltage the amplifier actually amplifies.  `vos` is
        ## the input voltage that NULLS the output, hence the sign.
        verr = _var(bd.V - vos + vcm / cmr                         # noqa
                    + (vsup - vsupnom) / psr, 'verr')              # noqa

        ## -- the gain node ----------------------------------------------
        ## Boyle's two relations, solved for the two element values.
        gm1 = _var(2.0 * sympy.pi * gbw * cc, 'gm1')               # noqa
        r2 = _var(aolc / gm1, 'r2')
        ## The slew-rate limit is a CURRENT limit on what charges `cc`,
        ## which is what it physically is in every compensated amplifier:
        ## the input stage's tail current is finite.  `tanh` because it
        ## is the smooth saturating function whose slope at the origin is
        ## exactly 1, so the small-signal transconductance is `gm1`
        ## unmodified and the gain relations above stay identities.  The
        ## relative error is `(gm1*verr/isl)^2/3`, which at the largest
        ## error voltage the amplifier ever sees in its linear region is
        ## below 1e-8.
        isl = _var(_maxc(sr * cc, 1e-30), 'islew')                 # noqa
        ig = _var(isl * sympy.tanh(gm1 * verr / isl), 'igain')
        b2 = Branch(Node('gain'), vee)
        ## `V(gain, vee)` is the amplifier's output SWING about the
        ## midpoint, not an absolute potential; the output stage adds the
        ## midpoint back.  That keeps the state the transient integrates
        ## centred on zero at quiescence.
        vg = _var(b2.V, 'vg')

        ## -- output stage ------------------------------------------------
        ## The available swing, and a smooth clamp to it.  `hypsmooth`
        ## rather than `tanh` because the clamp must not touch the gain
        ## in the linear region: `vs*tanh(v/vs)` is 8% low at half the
        ## swing, where this is `eps^2/(vs - v)` low -- 2e-6 V for a
        ## 15 V swing.  The smoothing width is written RELATIVE to the
        ## swing, since an absolute one would encode an assumption about
        ## a supply we do not control.
        vsw = _var(_maxc(vmid - vdrop, 1e-3), 'vswing')            # noqa
        esat = _var(1e-4 * vsw, 'esat')
        vclamp = _var(vsw - _hypsmooth(vsw - vg, esat)
                      + _hypsmooth(-vsw - vg, esat), 'vclamp')
        bo = Branch(Node('oint'), vee)
        bout = Branch(bo.plus, out)
        ## Output resistance AND short-circuit current from ONE
        ## expression: the slope at the origin is `1/rout` and the
        ## asymptote is `isc`.  Writing them as a resistor plus a
        ## separate clamp would need two elements and a corner between
        ## them; this is C-infinity and has no corner at all.
        routc = _var(_maxc(rout, 1e-6), 'routc')                   # noqa
        iscc = _var(_maxc(isc, 1e-12), 'iscc')                     # noqa

        return (
            Contribution(bd.I, bd.V / rin),                        # noqa
            Contribution(bp.I, (bp.V - vmid) / ricm + ib           # noqa
                         + 0.5 * ios),                             # noqa
            Contribution(bn.I, (bn.V - vmid) / ricm + ib           # noqa
                         - 0.5 * ios),                             # noqa
            Contribution(b2.I, vg / r2 + ddt(cc * vg) - ig),       # noqa
            Contribution(bo.V, vmid + vclamp),
            Contribution(bout.I, iscc * sympy.tanh(bout.V
                                                   / (routc * iscc))),
        )


## ======================================================================
## Fifth batch (2026-08-26): the VCO with real phase noise, the optical
## devices, the PLL's small behavioural blocks, the MESFET/HEMT family
## and MOS level 3.  Every class below is written with `params_as`: the
## parameters arrive as ``p`` and are read as attributes, so a helper
## shared by two classes takes ``p`` as an ordinary argument and a
## SPICE keyword-named parameter (``as``, ``lambda``) is declared under
## its own name and read as ``p['as']``.
## ======================================================================

_TWO_PI = 2 * sympy.pi


class VcoHdl(Behavioural):
    """A voltage-controlled oscillator whose noise is injected into the
    FREQUENCY, so that jitter accumulates as a random walk.  Terminals
    ``(cp, cn, outp, outn, ph)``::

        f          = f0 + kvco*V(cp,cn) + n_f(t)
        V(ph,outn) = idtmod(f, 0, modulus, 0)        # phase, in CYCLES
        V(outp,outn) = va*sin(2*pi*V(ph,outn))

    ``n_f`` is a white frequency noise of PSD ``sf`` (Hz^2/Hz) plus a
    flicker frequency noise ``kff/f`` (Hz^2).  Because it enters the
    integrand, the PHASE noise it produces is ``S_phi = sf/f_m^2``
    (``L(f_m)`` falling at -20 dB/decade) and ``kff/f_m^3`` -- Leeson's
    two slopes, from the integral itself and not from a shaped source.
    A noise source added to the output would be flat, and a source
    added to the phase would be flat in phase; neither accumulates.
    The roadmap (#6) calls this "a capability the repo's own survey
    found in no other simulator", and it is one line of DSL once the
    noise is where it belongs.

    **Where the noise actually lives, and why there is an internal
    node.**  ``idtmod(expr)`` refuses a noise term in ``expr``
    (`generate_code`: "noise inside idt/idtmod is not supported"), and
    a V-contribution refuses one too.  So the frequency is put on an
    internal node ``fn`` by a 1 S I-contribution -- ``I(fn,outn) <+
    V(fn,outn) - f`` -- and the noise is contributed as a CURRENT on
    that branch, which the 1 S turns into a frequency in Hz.  The phase
    integrates ``V(fn,outn)``.  It costs a node and a paragraph, and is
    recorded as friction (roadmap sec. 12, fifth batch).

    ``modulus`` is in cycles: the phase node wraps at ``modulus`` so a
    ``DividerHdl`` reading it with ``n = modulus`` sees a continuous
    ``sin(2*pi*V(ph)/n)``.  The output is periodic in one cycle
    whatever the modulus is.
    """
    params_as = 'p'
    instparams = [
        Parameter(name='f0', desc='Free-running frequency', unit='Hz',
                  default=1e6),
        Parameter(name='kvco', desc='Tuning gain', unit='Hz/V',
                  default=1e5),
        Parameter(name='va', desc='Output amplitude', unit='V', default=1.0),
        Parameter(name='modulus', desc='Phase wrap, in cycles', unit='',
                  default=1.0),
        Parameter(name='sf', desc='White frequency-noise PSD',
                  unit='Hz^2/Hz', default=0.0),
        Parameter(name='kff', desc='Flicker frequency-noise coefficient '
                  '(PSD kff/f)', unit='Hz^2', default=0.0),
    ]

    @staticmethod
    def analog(p, cp, cn, outp, outn, ph):
        bc, bo, bph = Branch(cp, cn), Branch(outp, outn), Branch(ph, outn)
        fn = Node('fn')
        bf = Branch(fn, outn)
        freq = _var(p.f0 + p.kvco * bc.V, 'freq')
        return (Contribution(bf.I, bf.V - freq),
                Contribution(bf.I, _white_noise(p.sf)),
                Contribution(bf.I, _flicker_noise(p.kff, 1)),
                Contribution(bph.V, idtmod(bf.V, 0.0, p.modulus, 0.0)),
                Contribution(bo.V, p.va * sympy.sin(_TWO_PI * bph.V)))


class DividerHdl(Behavioural):
    """A continuous ``/N`` frequency divider in the PHASE domain, and the
    second production user of `Cross`.  Terminals
    ``(inp, inn, outp, outn)``::

        s            = sin(2*pi*V(inp,inn)/n)
        V(outp,outn) = vol + (voh - vol)/2*(1 + tanh(gain*s))
        @cross(s, 0)

    The input is a PHASE in cycles -- `VcoHdl`'s ``ph`` terminal with
    ``modulus = n`` (or any multiple of ``n``), so that the wrap lands
    where ``s`` is already zero and the divided signal is continuous.
    Kundert's PLL modelling writes a divider exactly so
    (``phase_out = phase_in/N``); what is added here is the hard-limited
    output and the `Cross` at each of its edges, which is what a
    phase-frequency detector downstream would count.

    There is no counter, no latched state and no memory: the roadmap's
    "continuous" divider is the one that stays continuous by living in
    the phase domain.  A divider that counts EDGES of a real waveform
    needs event-latched integer state, which the DSL does not have
    (roadmap sec. 5, "What the behavioural half needs").
    """
    params_as = 'p'
    instparams = [
        Parameter(name='n', desc='Division ratio', unit='', default=4.0),
        Parameter(name='gain', desc='Output limiter gain', unit='',
                  default=50.0),
        Parameter(name='voh', desc='Output high', unit='V', default=1.0),
        Parameter(name='vol', desc='Output low', unit='V', default=-1.0),
    ]

    @staticmethod
    def analog(p, inp, inn, outp, outn):
        bi, bo = Branch(inp, inn), Branch(outp, outn)
        s = _var(sympy.sin(_TWO_PI * bi.V / p.n), 's')
        return (Contribution(bo.V, p.vol + (p.voh - p.vol) / 2
                             * (1 + sympy.tanh(p.gain * s))),
                Cross(s, 0))


class MixerHdl(Behavioural):
    """An ideal multiplying mixer.  Terminals
    ``(rfp, rfn, lop, lon, ifp, ifn)``::

        V(ifp,ifn) = k * V(rfp,rfn) * V(lop,lon)

    Two sines in give ``k*A*B/2`` at the difference and the sum
    frequency and nothing else -- the trigonometric identity is the
    whole model, and it is what the test asserts by projection.
    """
    params_as = 'p'
    instparams = [Parameter(name='k', desc='Conversion constant',
                            unit='1/V', default=1.0)]

    @staticmethod
    def analog(p, rfp, rfn, lop, lon, ifp, ifn):
        brf, blo, bif = Branch(rfp, rfn), Branch(lop, lon), Branch(ifp, ifn)
        return Contribution(bif.V, p.k * brf.V * blo.V)


class ChargePumpHdl(Behavioural):
    """A PLL charge pump.  Terminals ``(up, dn, outp, outn)``::

        I(outn -> outp) = iup*sw(V(up,outn)) - idn*sw(V(dn,outn))
        sw(v)           = (1 + tanh((v - vth)/vsw))/2

    UP sources ``iup`` INTO ``outp``; DN sinks ``idn`` from it; both on
    is the difference (a real pump's mismatch is ``iup != idn``), both
    off is zero to the last bit of ``tanh``.  ``vsw`` is the switching
    width of the tanh; a hard step would leave the controller no slope
    to find the edge with.
    """
    params_as = 'p'
    instparams = [
        Parameter(name='iup', desc='UP current', unit='A', default=100e-6),
        Parameter(name='idn', desc='DN current', unit='A', default=100e-6),
        Parameter(name='vth', desc='Logic threshold', unit='V',
                  default=0.5),
        Parameter(name='vsw', desc='Switching width', unit='V',
                  default=0.05),
    ]

    @staticmethod
    def analog(p, up, dn, outp, outn):
        bup, bdn, bo = Branch(up, outn), Branch(dn, outn), Branch(outp, outn)
        su = _var((1 + sympy.tanh((bup.V - p.vth) / p.vsw)) / 2, 'su')
        sd = _var((1 + sympy.tanh((bdn.V - p.vth) / p.vsw)) / 2, 'sd')
        ## `I(outp,outn) <+ x` is a current INTO the element at outp, so
        ## sourcing into the node is the negative.
        return Contribution(bo.I, -p.iup * su + p.idn * sd)


## ----------------------------------------------------------------------
## Optical devices -- roadmap item 13.  The optical port is a pair of
## terminals whose potential is an optical POWER in watts; a photodiode
## reads it and an LED drives it, and MNA does not care that the "volt"
## is a watt.  Both are the SPICE diode of `_spice_diode` with one
## extra statement.
## ----------------------------------------------------------------------

def _optical_params():
    return [
        Parameter(name='resp', desc='Responsivity', unit='A/W',
                  default=0.5),
        ## SPICE's "no shunt" is a large resistance, and 1e30 ohm is the
        ## arithmetic spelling of it: a 1e-30 S conductance is below the
        ## solver's abstol on any node it touches.
        Parameter(name='rsh', desc='Junction shunt resistance',
                  unit='ohm', default=1e30),
    ]


class PhotodiodeHdl(Behavioural):
    """A photodiode / solar cell: the SPICE diode with a photocurrent
    across the junction and a shunt resistance.  Terminals
    ``(a, c, optp, optn)``::

        I(c -> a, across the junction) = resp * V(optp,optn)
        I(junction)                   += V(junction)/rsh

    The photocurrent flows from cathode to anode inside the device --
    a lit photodiode sources current out of its anode -- and it sits
    BEHIND the series resistance, in parallel with the junction, which
    is the single-diode solar-cell circuit (Iph, D, Rsh, Rs).  With
    ``rs = 0`` and ``rsh = inf`` its I-V is the closed form
    ``I = Iph - IS*(exp(V/(n*Vt)) - 1)``, so ``Isc = Iph`` and
    ``Voc = n*Vt*ln(1 + Iph/IS)``, which the tests assert; with both
    parasitics the test solves the implicit form independently.

    A solar cell is this class with a cell's card (``IS`` of nA, ``n``
    around 1.3, a small ``rs``, a finite ``rsh``, ``resp`` scaled so
    that ``V(opt) = 1`` is one sun).  The parameters are SPICE's plus
    ``resp`` and ``rsh``; the full SPICE junction -- charge, breakdown,
    temperature, noise -- comes from `_spice_diode` unchanged.
    """
    params_as = 'p'
    instparams = _spice_diode_params() + _optical_params()

    @staticmethod
    def analog(p, a, c, optp, optn):
        bopt = Branch(optp, optn)
        iph = _var(p.resp * bopt.V, 'iph')
        stmts, _pd, _i = _spice_diode(
            p, a, c, TEMP, junction=lambda bd: -iph + bd.V / p.rsh)
        return stmts


class LedHdl(Behavioural):
    """A light-emitting diode: the SPICE diode with an optical output.
    Terminals ``(a, c, optp, optn)``::

        V(optp,optn) = eta * max(I_junction - ith, 0)

    ``eta`` is the slope efficiency in W/A and ``ith`` a threshold
    current (0 for an LED; a laser diode has one).  The optical branch
    is V-contributed, so it carries a branch-current unknown that a
    photodiode's optical port -- which only reads the potential --
    leaves at zero: the "current" through an optical node is nothing,
    and MNA is happy to solve for it.  `_spice_diode`'s junction
    current is the one reported to the port, so the light follows the
    forward current and not the displacement or breakdown current.
    """
    params_as = 'p'
    instparams = _spice_diode_params() + [
        Parameter(name='eta', desc='Slope efficiency', unit='W/A',
                  default=0.1),
        Parameter(name='ith', desc='Threshold current', unit='A',
                  default=0.0),
    ]

    @staticmethod
    def analog(p, a, c, optp, optn):
        bopt = Branch(optp, optn)
        stmts, _pd, idio = _spice_diode(p, a, c, TEMP)
        return stmts + (Contribution(bopt.V, p.eta * _maxc(idio - p.ith,
                                                            0.0)),)


## ----------------------------------------------------------------------
## MESFET / HEMT -- roadmap item 15.  Three tanh-shaped channels.
## ----------------------------------------------------------------------

def _mesfet_params():
    """SPICE's MESFET card (`mesdefs.h`), in SPICE's names -- including
    ``lambda``, declared as itself and read as ``p['lambda']``."""
    return [
        Parameter(name='vto', desc='Pinch-off voltage', unit='V',
                  default=-2.0),
        Parameter(name='beta', desc='Transconductance parameter',
                  unit='A/V^2', default=2.5e-3),
        Parameter(name='alpha', desc='Saturation voltage parameter',
                  unit='1/V', default=2.0),
        Parameter(name='lambda', desc='Channel-length modulation',
                  unit='1/V', default=0.0),
        Parameter(name='b', desc='Doping tail extending parameter',
                  unit='1/V', default=0.3),
        Parameter(name='rd', desc='Drain ohmic resistance', unit='ohm',
                  default=0.0),
        Parameter(name='rs', desc='Source ohmic resistance', unit='ohm',
                  default=0.0),
        Parameter(name='cgs', desc='Zero-bias G-S junction capacitance',
                  unit='F', default=0.0),
        Parameter(name='cgd', desc='Zero-bias G-D junction capacitance',
                  unit='F', default=0.0),
        Parameter(name='pb', desc='Gate junction potential', unit='V',
                  default=1.0),
        Parameter(name='IS', desc='Gate junction saturation current',
                  unit='A', default=1e-14),
        Parameter(name='fc', desc='Forward-bias depletion coefficient',
                  unit='', default=0.5),
        Parameter(name='kf', desc='Flicker-noise coefficient', unit='',
                  default=0.0),
        Parameter(name='af', desc='Flicker-noise exponent', unit='',
                  default=1.0),
        Parameter(name='area', desc='Area factor', unit='', default=1.0),
    ]


class MesfetStatzHdl(Behavioural):
    """The SPICE MESFET (Statz et al. 1987), as ngspice's `mesload.c`
    has it.  Terminals ``(d, g, s)``::

        Ids = beta*(1 + lambda*vds) * vgst^2/(1 + b*vgst)
              * [1 - (1 - alpha*vds/3)^3]     for 0 < vds < 3/alpha
              * 1                              for vds >= 3/alpha

    with ``vgst = max(vgs - vto, 0)``; for ``vds < 0`` the source and
    drain exchange roles (``vgd`` in place of ``vgs``, ``-vds``), which
    is written as a `Piecewise` on ``vds`` with the two arms meeting
    with equal value and slope at ``vds = 0``.  The gate is two
    Schottky junctions, ``IS*(exp(v/Vt) - 1)`` each, with SPICE's
    depletion charge (``cgs``/``cgd``/``pb``/``fc``), series ``rd`` /
    ``rs`` on collapsible internal nodes, channel thermal noise
    ``(8/3)*k*T*gm`` and flicker noise.

    Reference: H. Statz, P. Newman, I. Smith, R. Pucel and H. Haus,
    "GaAs FET device and circuit simulation in SPICE", *IEEE Trans.
    Electron Devices* **34**, 160 (1987); the ``b`` doping-tail term
    and the cubic saturation are theirs.  `mes1.va` (vacask's
    transcription of `mesload.c`) was the line-by-line reference.

    **The limiter declaration is NOT SPICE's, and the DSL is why.**
    `mesload.c` runs ``pnjlim`` AND ``fetlim`` on the same ``vgs``, in
    sequence; the DSL carries one law per branch (``branch is limited
    twice`` is a compile-time refusal).  The gate junctions are the
    exponentials, so they get ``pnjlim``, ``vds`` gets ``limvds``, and
    the three close a triangle, so they are one `limit_together` group.
    ``fetlim`` on the gate is what `MesfetCurticeHdl` and
    `HemtAngelovHdl` declare -- they have no gate junction.
    """
    params_as = 'p'
    instparams = _mesfet_params()

    @staticmethod
    def analog(p, d, g, s):
        di, si = Node('di'), Node('si')
        brd, brs = Branch(d, di), Branch(s, si)
        bgs, bgd, bds = Branch(g, si), Branch(g, di), Branch(di, si)
        vtT = _var(_vt(TEMP), 'vtT')
        isT = _var(p.area * p.IS, 'isT')
        vgs, vgd, vds = limit_together(limit_pnj(bgs.V, isT, vtT),
                                       limit_pnj(bgd.V, isT, vtT),
                                       limit_vds(bds.V))
        vgs, vgd, vds = _var(vgs, 'vgs'), _var(vgd, 'vgd'), _var(vds, 'vds')
        beta = _var(p.area * p.beta, 'betaa')
        lam = p['lambda']

        def channel(vg, vd, tag):
            vgst = _var(_maxc(vg - p.vto, 0.0), 'vgst' + tag)
            betap = _var(beta * (1 + lam * vd), 'betap' + tag)
            core = _var(betap * vgst * vgst / (1 + p.b * vgst), 'core' + tag)
            ## `1 - (1 - alpha*vd/3)^3` below `3/alpha`, 1 above -- and
            ## the cubic IS 1 at `vd = 3/alpha` with zero slope, so the
            ## `minc` is C1 rather than a clamp.
            afact = _var(1 - p.alpha * _minc(vd, 3.0 / p.alpha) / 3,
                         'afact' + tag)
            lfact = _var(1 - afact ** 3, 'lfact' + tag)
            return _var(core * lfact, 'ids' + tag)

        ifwd = channel(vgs, vds, 'f')
        irev = channel(vgd, -vds, 'r')
        ids = _var(sympy.Piecewise((ifwd, vds >= 0.0), (-irev, True)), 'ids')
        igs = _var(isT * (_expl(vgs / vtT) - 1), 'igs')
        igd = _var(isT * (_expl(vgd / vtT) - 1), 'igd')
        qgs = _var(_pn_depletion_charge(vgs, p.area * p.cgs, p.pb, 0.5,
                                        p.fc, 'gs'), 'qgs')
        qgd = _var(_pn_depletion_charge(vgd, p.area * p.cgd, p.pb, 0.5,
                                        p.fc, 'gd'), 'qgd')
        ## Channel noise: `mesnoise.c` uses `(8/3)*k*T*gm`.  For the
        ## square-law core `gm = 2*Ids/vgst` up to the doping-tail
        ## factor, which is what is used here; the floor keeps a
        ## cut-off device (vgst = 0, Ids = 0) at 0/0 = 0.
        gmn = _var(2 * _safe_abs(ids) / _maxc(
            sympy.Piecewise((vgs - p.vto, vds >= 0.0), (vgd - p.vto, True)),
            1e-3), 'gmn')
        stmts = (
            Contribution(bds.I, ids),
            Contribution(bgs.I, igs + ddt(qgs)),
            Contribution(bgd.I, igd + ddt(qgd)),
            Contribution(brd.I, brd.V * p.area / p.rd),
            Collapse(brd, p.rd <= 0.0),
            Contribution(brs.I, brs.V * p.area / p.rs),
            Collapse(brs, p.rs <= 0.0),
            Contribution(bds.I, _white_noise(8.0 / 3.0 * _KB * TEMP * gmn)),
            Contribution(bds.I, _flicker_noise(p.kf * _safe_abs(ids) ** p.af,
                                               1)),
            Contribution(brd.I, _white_noise(4 * _KB * TEMP * p.area / p.rd)),
            Contribution(brs.I, _white_noise(4 * _KB * TEMP * p.area / p.rs)),
        )
        return stmts


class MesfetCurticeHdl(Behavioural):
    """Curtice's quadratic MESFET (1980).  Terminals ``(d, g, s)``::

        Ids = beta * max(vgs - vto, 0)^2 * (1 + lambda*vds) * tanh(alpha*vds)

    W. R. Curtice, "A MESFET model for use in the design of GaAs
    integrated circuits", *IEEE Trans. Microwave Theory Tech.* **28**,
    448 (1980), eq. (1).  Intrinsic channel only -- no gate junction,
    no charge -- which is the equation as published and the one whose
    ``tanh`` gives the `Statz` model its saturation limit at
    ``alpha -> infinity``.  The ``tanh`` is odd, so the equation is
    antisymmetric in ``vds`` at fixed ``vgs`` (it does not exchange
    source and drain; Curtice's paper does not either).

    ``limit_fet`` on the gate and ``limit_vds`` on the drain, as one
    `limit_together` group: the channel is a polynomial times a bounded
    function and cannot overflow, so what the limiters buy here is
    measured rather than assumed (`test_elements_hdl_library5.py`).
    """
    params_as = 'p'
    instparams = [
        Parameter(name='vto', desc='Pinch-off voltage', unit='V',
                  default=-2.0),
        Parameter(name='beta', desc='Transconductance parameter',
                  unit='A/V^2', default=2.5e-3),
        Parameter(name='alpha', desc='Saturation voltage parameter',
                  unit='1/V', default=2.0),
        Parameter(name='lambda', desc='Channel-length modulation',
                  unit='1/V', default=0.0),
    ]

    @staticmethod
    def analog(p, d, g, s):
        bgs, bds = Branch(g, s), Branch(d, s)
        vgs, vds = limit_together(limit_fet(bgs.V, p.vto), limit_vds(bds.V))
        vgs, vds = _var(vgs, 'vgs'), _var(vds, 'vds')
        vgst = _var(_maxc(vgs - p.vto, 0.0), 'vgst')
        ids = _var(p.beta * vgst * vgst * (1 + p['lambda'] * vds)
                   * sympy.tanh(p.alpha * vds), 'ids')
        return Contribution(bds.I, ids)


class HemtAngelovHdl(Behavioural):
    """Angelov's HEMT/MESFET model (1992).  Terminals ``(d, g, s)``::

        Ids   = ipk * (1 + tanh(psi)) * (1 + lambda*vds) * tanh(alpha*vds)
        psi   = p1*(vgs - vpk) + p2*(vgs - vpk)^2 + p3*(vgs - vpk)^3
        alpha = alphar + alphas*(1 + tanh(psi))

    I. Angelov, H. Zirath and N. Rorsman, "A new empirical nonlinear
    model for HEMT and MESFET devices", *IEEE Trans. Microwave Theory
    Tech.* **40**, 2258 (1992), eqs. (1)-(3).  ``ipk`` and ``vpk`` are
    the current and gate voltage at PEAK transconductance: with
    ``p2 = p3 = 0``, ``gm = ipk*p1*sech^2(psi)`` is largest at
    ``vgs = vpk``, where ``Ids = ipk`` (saturated) -- both identities
    of the parameterisation, and both asserted.  Intrinsic channel
    only: Angelov's ``Cgs``/``Cgd`` are capacitance fits, not the
    gradient of a charge, and the DSL contributes charge (the same
    reason `MosLevel1Hdl` omits Meyer).

    ``limit_fet(vgs, vpk)`` -- the law measures excursions against the
    peak-gm voltage, the nearest thing this model has to a threshold --
    and ``limit_vds``, grouped.
    """
    params_as = 'p'
    instparams = [
        Parameter(name='ipk', desc='Drain current at peak gm', unit='A',
                  default=50e-3),
        Parameter(name='vpk', desc='Gate voltage at peak gm', unit='V',
                  default=-0.2),
        Parameter(name='p1', desc='psi first-order coefficient',
                  unit='1/V', default=1.5),
        Parameter(name='p2', desc='psi second-order coefficient',
                  unit='1/V^2', default=0.0),
        Parameter(name='p3', desc='psi third-order coefficient',
                  unit='1/V^3', default=0.0),
        Parameter(name='alphar', desc='Saturation parameter, residual',
                  unit='1/V', default=1.0),
        Parameter(name='alphas', desc='Saturation parameter, gate-'
                  'dependent part', unit='1/V', default=0.5),
        Parameter(name='lambda', desc='Channel-length modulation',
                  unit='1/V', default=0.02),
    ]

    @staticmethod
    def analog(p, d, g, s):
        bgs, bds = Branch(g, s), Branch(d, s)
        vgs, vds = limit_together(limit_fet(bgs.V, p.vpk), limit_vds(bds.V))
        vgs, vds = _var(vgs, 'vgs'), _var(vds, 'vds')
        dv = _var(vgs - p.vpk, 'dv')
        psi = _var(p.p1 * dv + p.p2 * dv ** 2 + p.p3 * dv ** 3, 'psi')
        tp = _var(1 + sympy.tanh(psi), 'tp')
        alpha = _var(p.alphar + p.alphas * tp, 'alpha')
        ids = _var(p.ipk * tp * (1 + p['lambda'] * vds)
                   * sympy.tanh(alpha * vds), 'ids')
        return Contribution(bds.I, ids)


## ----------------------------------------------------------------------
## SPICE MOS level 3 -- roadmap item 10's second half
## ----------------------------------------------------------------------

def _mos3_params():
    """The SPICE level-3 card, in SPICE's own names -- ``as`` included,
    which the namespace lets a model declare outright."""
    return [
        Parameter(name='vto', desc='Zero-bias threshold voltage', unit='V',
                  default=0.0),
        Parameter(name='kp', desc='Transconductance parameter (0: from '
                  'u0 and tox)', unit='A/V^2', default=0.0),
        Parameter(name='u0', desc='Surface mobility', unit='cm^2/V/s',
                  default=600.0),
        Parameter(name='gamma', desc='Body effect factor (unset: from '
                  'nsub and tox)', unit='V^0.5', default=0.0),
        Parameter(name='phi', desc='Surface inversion potential (unset: '
                  'from nsub)', unit='V', default=0.6),
        Parameter(name='tox', desc='Oxide thickness', unit='m',
                  default=1e-7),
        Parameter(name='nsub', desc='Substrate doping', unit='cm^-3',
                  default=0.0),
        Parameter(name='xj', desc='Junction depth', unit='m', default=0.0),
        Parameter(name='nfs', desc='Fast surface state density',
                  unit='cm^-2', default=0.0),
        Parameter(name='eta', desc='Static feedback (DIBL) coefficient',
                  unit='', default=0.0),
        Parameter(name='delta', desc='Width effect on threshold', unit='',
                  default=0.0),
        Parameter(name='theta', desc='Mobility degradation', unit='1/V',
                  default=0.0),
        Parameter(name='vmax', desc='Maximum drift velocity (0: none)',
                  unit='m/s', default=0.0),
        Parameter(name='kappa', desc='Saturation field factor', unit='',
                  default=0.2),
        Parameter(name='w', desc='Channel width', unit='m', default=1e-4),
        Parameter(name='l', desc='Drawn channel length', unit='m',
                  default=1e-4),
        Parameter(name='ld', desc='Lateral diffusion', unit='m',
                  default=0.0),
        Parameter(name='wd', desc='Width reduction', unit='m', default=0.0),
        Parameter(name='cgso', desc='Gate-source overlap capacitance per '
                  'metre', unit='F/m', default=0.0),
        Parameter(name='cgdo', desc='Gate-drain overlap capacitance per '
                  'metre', unit='F/m', default=0.0),
        Parameter(name='cgbo', desc='Gate-bulk overlap capacitance per '
                  'metre', unit='F/m', default=0.0),
        Parameter(name='cbd', desc='Zero-bias bulk-drain capacitance '
                  '(0: cj*ad)', unit='F', default=0.0),
        Parameter(name='cbs', desc='Zero-bias bulk-source capacitance '
                  '(0: cj*as)', unit='F', default=0.0),
        Parameter(name='IS', desc='Bulk junction saturation current',
                  unit='A', default=1e-14),
        Parameter(name='pb', desc='Bulk junction potential', unit='V',
                  default=0.8),
        Parameter(name='cj', desc='Bulk bottom capacitance per area',
                  unit='F/m^2', default=0.0),
        Parameter(name='cjsw', desc='Bulk sidewall capacitance per metre',
                  unit='F/m', default=0.0),
        Parameter(name='mj', desc='Bottom grading coefficient', unit='',
                  default=0.5),
        Parameter(name='mjsw', desc='Sidewall grading coefficient',
                  unit='', default=0.33),
        Parameter(name='fc', desc='Forward-bias depletion coefficient',
                  unit='', default=0.5),
        Parameter(name='js', desc='Bulk junction saturation current per '
                  'area', unit='A/m^2', default=0.0),
        Parameter(name='rd', desc='Drain ohmic resistance', unit='ohm',
                  default=0.0),
        Parameter(name='rs', desc='Source ohmic resistance', unit='ohm',
                  default=0.0),
        Parameter(name='rsh', desc='Sheet resistance', unit='ohm/sq',
                  default=0.0),
        Parameter(name='ad', desc='Drain diffusion area', unit='m^2',
                  default=0.0),
        Parameter(name='as', desc='Source diffusion area', unit='m^2',
                  default=0.0),
        Parameter(name='pd', desc='Drain diffusion perimeter', unit='m',
                  default=0.0),
        Parameter(name='ps', desc='Source diffusion perimeter', unit='m',
                  default=0.0),
        Parameter(name='nrd', desc='Drain squares', unit='', default=1.0),
        Parameter(name='nrs', desc='Source squares', unit='', default=1.0),
        Parameter(name='kf', desc='Flicker-noise coefficient', unit='',
                  default=0.0),
        Parameter(name='af', desc='Flicker-noise exponent', unit='',
                  default=1.0),
        Parameter(name='tnom', desc='Parameter measurement temperature',
                  unit='K', default=300.15),
    ]


def _mos3_channel(p, vgs, vds, vbs, c, tag):
    """One ARM of the level-3 channel: ``vds >= 0``, in the n-channel
    convention.  ``c`` carries the card-derived constants (a dict);
    returns ``(cdrain, von, vdsat)``, every intermediate named with
    ``tag`` so the two arms coexist in one chain.

    Transcribed from `mos3load.c` through `mos3.va` (vacask's
    Verilog-A distillation of ngspice-44), in its order: the
    forward-biased bulk continuation, ``fshort`` (Dang's charge-sharing
    with SPICE's three coefficients), the body factor, ``eta``,
    ``delta``, the ``nfs`` subthreshold slope, ``theta``, ``vmax``
    saturation, the ``kappa`` channel-length modulation in its three
    cases and the ``delxl > Leff/2`` continuation, the subthreshold
    exponential.  Every ``if`` on the solution is a `Piecewise`;
    every ``if`` on a parameter is a `Piecewise` on the parameter with
    BOTH arms kept finite, because both are evaluated.
    """
    phiT, leff, weff, cox, vtT = c['phiT'], c['leff'], c['weff'], c['cox'], \
        c['vtT']
    gam, vbiT, kpT, u0T = c['gam'], c['vbiT'], c['kpT'], c['u0T']

    def v(expr, name):
        return _var(expr, name + tag)

    ## -- the bulk continuation (SPICE's, also Level 1's) --------------
    ## The bulk continuation, clamped by `hdl.select` from each arm's
    ## own condition rather than by a chosen floor; see `sarg` in
    ## `_spice_mos1` for the argument and for what it replaces.
    sqphbs = v(_select(
        (sympy.sqrt(phiT - vbs), vbs <= 0.0),
        (sympy.sqrt(phiT) / (1.0 + vbs / (2.0 * phiT)), True)),
        'sqphbs')
    phibs = v(sqphbs * sqphbs, 'phibs')

    ## -- fshort: short-channel charge sharing ---------------------------
    ## `xd` (SPICE's coeffDepLayWidth) is sqrt(2*eps/(q*nsub)); zero
    ## without nsub, and then fshort is 1.  Every divisor floored so
    ## the discarded arm stays finite at `xj = 0`.
    xjf = v(_maxc(p.xj, 1e-12), 'xjf')
    wps = v(c['xd'] * sqphbs, 'wps')
    wponxj = v(wps / xjf, 'wponxj')
    wconxj = v(0.0631353 + 0.8013292 * wponxj - 0.01110777 * wponxj ** 2,
               'wconxj')
    arga = v(wconxj + p.ld / xjf, 'argafs')
    argc = v(wponxj / (1.0 + wponxj), 'argcfs')
    argb = v(sympy.sqrt(_maxc(1.0 - argc * argc, 1e-12)), 'argbfs')
    fshort = v(sympy.Piecewise(
        (1.0 - p.xj / leff * (arga * argb - p.ld / xjf),
         sympy.And(p.xj > 0.0, p.nsub > 0.0)),
        (1.0, True)), 'fshort')

    ## -- body factor, threshold ---------------------------------------
    gammas = v(gam * fshort, 'gammas')
    fbody = v(0.5 * gammas / (2.0 * sqphbs) + c['narrow'] / weff, 'fbody')
    qbonco = v(gammas * sqphbs + c['narrow'] * phibs / weff, 'qbonco')
    vth = v(vbiT - c['etal'] * vds + qbonco, 'vth')
    ## -- subthreshold slope factor (nfs) ------------------------------
    xn = v(1.0 + c['csonco'] + qbonco / (2.0 * phibs), 'xn')
    von = v(sympy.Piecewise((vth + vtT * xn, p.nfs > 0.0), (vth, True)),
            'von')
    vgsx = v(_maxc(vgs, von), 'vgsx')
    ## -- mobility degradation and the saturation voltage --------------
    onfg = v(1.0 + p.theta * (vgsx - vth), 'onfg')
    fgate = v(1.0 / onfg, 'fgate')
    us = v(u0T * 1e-4 * fgate, 'us')
    vdsc = v(leff * _maxc(p.vmax, 1e-30) / us, 'vdsc')
    arg0 = v((vgsx - vth) / (1.0 + fbody), 'arg0')
    vdsat = v(sympy.Piecewise(
        (arg0 + vdsc - sympy.sqrt(arg0 * arg0 + vdsc * vdsc), p.vmax > 0.0),
        (arg0, True)), 'vdsat')
    ## -- the current ---------------------------------------------------
    vdsx = v(_minc(vds, vdsat), 'vdsx')
    cdo = v(vgsx - vth - 0.5 * (1.0 + fbody) * vdsx, 'cdo')
    beta = v(kpT * weff / leff * fgate, 'beta')
    fdrain = v(sympy.Piecewise((1.0 / (1.0 + vdsx / vdsc), p.vmax > 0.0),
                               (1.0, True)), 'fdrain')
    cd0 = v(beta * cdo * vdsx * fdrain, 'cd0')
    ## -- channel-length modulation (kappa) ------------------------------
    ## Three cases, as `mos3load.c`: below vdsat without vmax the
    ## quartic ramp; above vdsat without vmax the plain root; above
    ## vdsat with vmax the `emax` form.  With vmax and below vdsat, or
    ## with no nsub (alpha = 0), there is none.  The ratio `vds/vdsat`
    ## is clamped at 1 -- it is only SELECTED below vdsat, where it is
    ## at most 1, and unclamped it is a 1e120 in the discarded arm.
    vdsatf = v(_maxc(vdsat, 1e-9), 'vdsatf')
    ratio = v(_minc(vds / vdsatf, 1.0), 'ratio')
    ka = v(p.kappa * c['alpha'], 'ka')
    ## `sqrt(ka*vdsat/8)` is differentiated as `sqrt(u)/(2u)`, which is
    ## 0/0 at `vdsat = 0` -- every cut-off device's reverse arm.  The
    ## floor is twenty decades below any selected value.
    dl_below = v(sympy.sqrt(_maxc(ka * vdsat, 1e-30) / 8.0) * ratio ** 4,
                 'dlbelow')
    ## Floored at 1e-30, NOT at 0: `sqrt(maxc(u, 0))` has derivative
    ## `inf * 0 = NaN` wherever `u < 0`, and this arm is evaluated in
    ## the reverse mode with `vds < 0` at every bias.  Selected, `u`
    ## is at least `ka*vdsat/8 > 0` and the floor is never reached.
    dl_novmax = v(sympy.sqrt(_maxc(ka * (vds - vdsat + vdsat / 8.0), 1e-30)),
                  'dlnovmax')
    cdsat = v(_maxc(cd0, 1e-30), 'cdsat')
    gdsat = v(_maxc(cdsat * (1.0 - fdrain) / vdsc, 1e-12), 'gdsat')
    emax = v(p.kappa * cdsat / (leff * gdsat), 'emax')
    ea = v(0.5 * emax * c['alpha'], 'ea')
    dl_vmax = v(sympy.sqrt(ea * ea + _maxc(ka * (vds - vdsat), 0.0) + 1e-60)
                - ea, 'dlvmax')
    delxl0 = v(sympy.Piecewise(
        (0.0, c['alpha'] <= 0.0),
        (sympy.Piecewise((0.0, p.vmax > 0.0), (dl_below, True)),
         vds <= vdsat),
        (dl_vmax, p.vmax > 0.0),
        (dl_novmax, True)), 'delxl0')
    ## The hand-written floor here was `_maxc(delxl0, 1e-30)` -- a
    ## constant, chosen to be below anything the arm is selected at.
    ## `select` derives the floor from the CONDITION instead: the arm
    ## is selected where `delxl0 > 0.5*leff`, so that is the floor, and
    ## `4*maxc(delxl0, 0.5*leff) >= 2*leff` needs no constant at all.
    ## The two agree bit for bit wherever the arm is selected (both
    ## return `delxl0` there, with derivative exactly 1) and differ
    ## only in the discarded arm, which is the point.
    delxl = v(_select(
        (leff - leff * leff / (4.0 * delxl0), delxl0 > 0.5 * leff),
        (delxl0, True)), 'delxl')
    xlfact = v(1.0 / (1.0 - delxl / leff), 'xlfact')
    ## -- subthreshold ----------------------------------------------------
    wfact = v(sympy.Piecewise(
        (_expl((vgs - von) / (vtT * _maxc(xn, 1.0))), vgs < von),
        (1.0, True)), 'wfact')
    cdrain = v(cd0 * xlfact * wfact, 'cdrain')
    return cdrain, von, vdsat


def _mos3_analog(T, nmos):
    """Build the ``analog()`` body of a SPICE level-3 MOSFET.

    Equations: `mos3load.c` / `mos3temp.c` (Berkeley SPICE3, ngspice-44
    through vacask's `mos3.va`), which is the reference the tests
    transcribe independently; the text is Massobrio & Antognetti,
    *Semiconductor Device Modeling with SPICE*, 2nd ed., sec. 3.4.
    The level-3 additions over level 1 are ``xj`` (short-channel
    threshold), ``delta`` (narrow-channel), ``eta`` (DIBL), ``nfs``
    (subthreshold), ``theta`` (mobility), ``vmax`` (velocity
    saturation) and ``kappa`` (channel-length modulation); the
    bulk-charge form of the current -- ``(1 + fbody)/2`` in the triode
    term and ``vdsat = (vgs - vth)/(1 + fbody)`` -- is level 3's own,
    and it coincides with level 1 only at ``gamma = 0``, which is the
    limit the tests use.

    **Source/drain exchange is a `Piecewise` on ``vds``**, unlike
    level 1's single antisymmetric expression: the level-3 channel has
    ``vdsat``, ``vmax`` and ``kappa`` in it and no identity makes it
    antisymmetric.  The two arms share every card-derived quantity and
    are built by `_mos3_channel` with a tag each; they meet at
    ``vds = 0`` with equal value (0) and equal slope
    (``beta*(vgs - vth)``).  It doubles the channel chain, and that is
    the compile-time price of a model written the way SPICE branches.

    Everything else -- both bulk junctions' current and depletion
    charge, the overlap capacitances, ``rd``/``rs`` on collapsible
    nodes, ``gamma``/``phi`` given-or-derived, the temperature path,
    the noise -- is level 1's, restated through ``p``.  Meyer is absent
    for level 1's reason.
    """
    def analog(p, d, g, s, b):
        di, si = Node('di'), Node('si')
        brd, brs = Branch(d, di), Branch(s, si)
        if nmos > 0:
            bgs, bds = Branch(g, si), Branch(di, si)
            bbs, bbd = Branch(b, si), Branch(b, di)
            bgd, bgb = Branch(g, di), Branch(g, b)
        else:
            bgs, bds = Branch(si, g), Branch(si, di)
            bbs, bbd = Branch(si, b), Branch(di, b)
            bgd, bgb = Branch(di, g), Branch(b, g)

        ## -- geometry and the derived card --------------------------
        leff = _var(_maxc(p.l - 2.0 * p.ld, 1e-9), 'leff')
        weff = _var(_maxc(p.w - 2.0 * p.wd, 1e-9), 'weff')
        cox = _var(_EPSOX / _maxc(p.tox, 1e-12), 'cox')
        vtn = _var(_KB * p.tnom / _QE, 'vtnom')
        nsm3 = _var(_maxc(p.nsub, 0.0) * 1e6, 'nsm3')
        gd = _var(sympy.sqrt(2.0 * _QE * _EPSSI * nsm3) / cox, 'gamd')
        phid = _var(2.0 * vtn * sympy.log(_maxc(nsm3 / (_NI_CM3 * 1e6),
                                                1.0)), 'phid')
        gg = _var(p.given('gamma'), 'ggam')
        gp = _var(p.given('phi'), 'gphi')
        gam = _var(gg * p.gamma + (1.0 - gg) * gd, 'gam')
        phz = _var(_maxc(gp * p.phi + (1.0 - gp) * phid, 0.1), 'phz')
        ## SPICE's `alpha = 2*eps_si/(q*nsub)` and `xd = sqrt(alpha)`:
        ## zero without `nsub`, which switches fshort and kappa off.
        ## `select`, not `Piecewise`: the arm is selected exactly where
        ## `nsm3 > 0`, so the clamp that keeps the OTHER arm's evaluation
        ## finite is derivable from the condition instead of hand-picked.
        ## `margin=1.0` reproduces the previous `_maxc(nsm3, 1.0)` bit for
        ## bit -- `nsm3` is a doping density near 1e21, so a 1.0 floor is
        ## far below anything selected.
        alpha = _var(_select((2.0 * _EPSSI / (_QE * nsm3), nsm3 > 0.0),
                             (0.0, True), margin=1.0), 'alpha')
        xd = _var(sympy.sqrt(alpha), 'xd')
        narrow = _var(p.delta * 0.5 * sympy.pi * _EPSSI / cox, 'narrow')
        etal = _var(p.eta * 8.15e-22 / (cox * leff ** 3), 'etal')
        csonco = _var(_QE * p.nfs * 1e4 / cox, 'csonco')
        kp0 = _var(sympy.Piecewise((p.kp, p.kp > 0.0),
                                   (p.u0 * cox * 1e-4, True)), 'kp0')

        ## -- temperature (`mos3temp.c`) ------------------------------
        vtT = _var(_vt(T), 'vtT')
        trat = _var(T / p.tnom, 'trat')
        ltr = _var(sympy.log(trat), 'ltrat')
        egT = _var(1.16 - 7.02e-4 * T ** 2 / (T + 1108.0), 'egT')
        egn = _var(1.16 - 7.02e-4 * p.tnom ** 2 / (p.tnom + 1108.0),
                   'egtnom')
        phiT = _var(_maxc(phz * trat - 3.0 * vtT * ltr
                          - egn * trat + egT, 0.1), 'phiT')
        pbT = _var(p.pb * trat - 3.0 * vtT * ltr - egn * trat + egT, 'pbT')
        ## `tVbi = vto - gamma*sqrt(phi) + (Eg(tnom) - Eg(T))/2
        ##        + (phi(T) - phi)/2`, then `vth(T) = tVbi + gamma*sqrt
        ## (phi(T))`.  The `(Eg(tnom) - Eg(T))/2` term is in mos1temp.c
        ## and mos3temp.c alike; `MosLevel1Hdl` above omits it (its
        ## docstring's `vbi(T)` line is incomplete).  Recorded, not
        ## changed there: it is a `T != tnom` effect of ~0.5*dEg only.
        vbiT = _var(p.vto - gam * sympy.sqrt(phz) + 0.5 * (egn - egT)
                    + 0.5 * (phiT - phz), 'vbiT')
        ratio4 = _var(_safe_pow(trat, -1.5, lo=1e-3), 'ratio4')
        kpT = _var(kp0 * ratio4, 'kpT')
        u0T = _var(p.u0 * ratio4, 'u0T')
        factlog = _var(-egT / vtT + egn / vtn, 'factlog')
        tshift = _var(4e-4 * (T - p.tnom) - (pbT / p.pb - 1.0), 'tshift')
        cfacb = _var(1.0 + p.mj * tshift, 'cfacb')
        cfacw = _var(1.0 + p.mjsw * tshift, 'cfacw')
        isbd = _var(_maxc(sympy.Piecewise((p.js * p.ad, p.js * p.ad > 0.0),
                                          (p.IS, True))
                          * _expl(factlog), 1e-30), 'isbd')
        isbs = _var(_maxc(sympy.Piecewise((p.js * p['as'],
                                           p.js * p['as'] > 0.0),
                                          (p.IS, True))
                          * _expl(factlog), 1e-30), 'isbs')

        ## -- the four probes, as one group (SPICE's set) ------------
        vbsl = _var(bbs.V, 'vbsl')
        sargl = _var(_select(
            (sympy.sqrt(phiT - vbsl), vbsl <= 0.0),
            (sympy.sqrt(phiT) / (1.0 + vbsl / (2.0 * phiT)),
             True)), 'sargl')
        vonl = _var(vbiT + gam * sargl, 'vonl')
        vgs, vds, vbs, vbd = limit_together(
            limit_fet(bgs.V, vonl), limit_vds(bds.V),
            limit_pnj(bbs.V, isbs, vtT), limit_pnj(bbd.V, isbd, vtT))
        vgs, vds = _var(vgs, 'vgs'), _var(vds, 'vds')
        vbs, vbd = _var(vbs, 'vbs'), _var(vbd, 'vbd')
        vgd = _var(vgs - vds, 'vgd')

        ## -- the channel, both modes --------------------------------
        cst = dict(phiT=phiT, leff=leff, weff=weff, cox=cox, vtT=vtT,
                   gam=gam, vbiT=vbiT, kpT=kpT, u0T=u0T, xd=xd,
                   narrow=narrow, etal=etal, csonco=csonco, alpha=alpha)
        idf, vonf, vdsatf = _mos3_channel(p, vgs, vds, vbs, cst, 'f')
        idr, _vonr, _vdsatr = _mos3_channel(p, vgd, -vds, vbd, cst, 'r')
        ids = _var(sympy.Piecewise((idf, vds >= 0.0), (-idr, True)), 'ids')

        ## -- the two bulk junctions ---------------------------------
        ibs = _var(isbs * (_expl(vbs / vtT) - 1.0), 'ibs')
        ibd = _var(isbd * (_expl(vbd / vtT) - 1.0), 'ibd')
        cbdb = _var(cfacb * sympy.Piecewise((p.cbd, p.cbd > 0.0),
                                            (p.cj * p.ad, True)), 'cbdb')
        cbsb = _var(cfacb * sympy.Piecewise((p.cbs, p.cbs > 0.0),
                                            (p.cj * p['as'], True)), 'cbsb')
        cbdw = _var(cfacw * p.cjsw * p.pd, 'cbdw')
        cbsw = _var(cfacw * p.cjsw * p.ps, 'cbsw')
        qbd = _var(_pn_depletion_charge(vbd, cbdb, pbT, p.mj, p.fc, 'bd')
                   + _pn_depletion_charge(vbd, cbdw, pbT, p.mjsw, p.fc,
                                          'bdw'), 'qbd')
        qbs = _var(_pn_depletion_charge(vbs, cbsb, pbT, p.mj, p.fc, 'bs')
                   + _pn_depletion_charge(vbs, cbsw, pbT, p.mjsw, p.fc,
                                          'bsw'), 'qbs')
        rdx = _var(sympy.Piecewise((p.rd, p.rd > 0.0), (p.rsh * p.nrd, True)),
                   'rdx')
        rsx = _var(sympy.Piecewise((p.rs, p.rs > 0.0), (p.rsh * p.nrs, True)),
                   'rsx')

        ## -- noise: (8/3)kT gm for a saturated channel, Nyquist for a
        ## linear one, through the same bulk-charge form as the current
        ## (`mos3noi.c` uses `(8/3)*k*T*(gm + gds + gmbs)`; here the
        ## Klaassen-Prins integral over the level-3 triode expression).
        vgtn = _var(_maxc(sympy.Piecewise((vgs, vds >= 0.0), (vgd, True))
                          - sympy.Piecewise((vonf, vds >= 0.0),
                                            (_vonr, True)), 0.0), 'vgtn')
        vdsn = _var(_minc(_safe_abs(vds), sympy.Piecewise(
            (vdsatf, vds >= 0.0), (_vdsatr, True))), 'vdsn')
        gn = _var(_safe_div(vgtn * vgtn - vgtn * vdsn + vdsn * vdsn / 3.0,
                            _maxc(vgtn - 0.5 * vdsn, 1e-9)) * kpT
                  * weff / leff, 'gn')

        stmts = (
            Contribution(bds.I, ids),
            Contribution(bbs.I, ibs + ddt(qbs)),
            Contribution(bbd.I, ibd + ddt(qbd)),
            Contribution(bgs.I, ddt(p.cgso * weff * bgs.V)),
            Contribution(bgd.I, ddt(p.cgdo * weff * bgd.V)),
            Contribution(bgb.I, ddt(p.cgbo * leff * bgb.V)),
            Contribution(brd.I, brd.V / rdx),
            Collapse(brd, sympy.And(p.rd <= 0.0, p.rsh * p.nrd <= 0.0)),
            Contribution(brs.I, brs.V / rsx),
            Collapse(brs, sympy.And(p.rs <= 0.0, p.rsh * p.nrs <= 0.0)),
            Contribution(bds.I, _white_noise(4.0 * _KB * T * gn)),
            Contribution(bds.I, _flicker_noise(
                p.kf * _safe_abs(ids) ** p.af / (cox * leff ** 2), 1)),
            Contribution(brd.I, _white_noise(4.0 * _KB * T / rdx)),
            Contribution(brs.I, _white_noise(4.0 * _KB * T / rsx)),
        )
        return stmts
    return analog


class MosLevel3Hdl(Behavioural):
    """SPICE level-3 n-channel MOSFET.  Terminals ``(d, g, s, b)``.

    Parameters and defaults are SPICE's, ``as`` and all -- the first
    model in the library to declare a keyword-named parameter under
    its own name.  See `_mos3_analog` for the equations, the reference
    and the source/drain `Piecewise`; `_mos3_channel` is the channel.
    Limiting is level 1's four-probe `limit_together` group.

    **Partial, and where.**  Shipped: the whole DC channel (``xj``,
    ``delta``, ``eta``, ``nfs``, ``theta``, ``vmax``, ``kappa``, both
    modes), both bulk junctions with depletion charge, overlap
    capacitance, ``rd``/``rs``, the temperature path, channel and
    flicker noise.  Absent: the Meyer intrinsic capacitances (level 1's
    reason), ``xl``/``xw`` geometry offsets and ``delvto``.
    """
    params_as = 'p'
    instparams = _mos3_params()

    analog = staticmethod(_mos3_analog(TEMP, +1))


class MosLevel3PmosHdl(Behavioural):
    """SPICE level-3 p-channel MOSFET.  Terminals ``(d, g, s, b)``;
    `MosLevel3Hdl` with every branch reversed, positive magnitudes on
    the card (see `MosLevel1PmosHdl`)."""
    params_as = 'p'
    instparams = _mos3_params()

    analog = staticmethod(_mos3_analog(TEMP, -1))
