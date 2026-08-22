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
                                   ddt, idt, idtmod, TIME)


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
