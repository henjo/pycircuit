"""The HDL element classes of the shooting tests.  Kept in a file of
their own: the HDL cache hashes the file an `analog` is defined in, so
an edit anywhere in a shared file would recompile them all.
"""
from pycircuit.circuit import *
from pycircuit.circuit.shooting import (PAC, algebraic_conditioning,
                                        topological_index)
import warnings
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution,
                                   Parameter as _HdlParameter, white_noise)
from pycircuit.post import Waveform, average
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal
import unittest
import pytest
import functools as _functools


class _SwitchHdl(Behavioural):
    """A behavioural switch mirroring a Verilog-A `pcswitch`, line for line:
    `g = goff + (gon - goff) * (1 + tanh((V(cp,cn) - vth)/vs))/2`,
    `I(p,n) <+ g V(p,n) + white_noise(4 kb T g)`.  `kb` and `temp` are
    parameters so nothing hides in a constant.  The noise is
    CYCLOSTATIONARY by construction."""
    terminals = ('p', 'n', 'cp', 'cn')
    instparams = [
        _HdlParameter(name='gon', desc='Closed conductance', unit='S',
                      default=1e-3),
        _HdlParameter(name='goff', desc='Open conductance', unit='S',
                      default=1e-9),
        _HdlParameter(name='vth', desc='Gate threshold', unit='V',
                      default=0.0),
        _HdlParameter(name='vs', desc='Softening', unit='V', default=50e-3),
        _HdlParameter(name='temp', desc='Noise temperature', unit='K',
                      default=300.0),
        _HdlParameter(name='kb', desc='Boltzmann constant', unit='J/K',
                      default=1.38e-23)]

    @staticmethod
    def analog(p, n, cp, cn):
        import sympy
        b = Branch(p, n)
        ctrl = Branch(cp, cn)
        s = (1 + sympy.tanh((ctrl.V - vth) / vs)) / 2          # noqa: F821
        g = goff + (gon - goff) * s                            # noqa: F821
        return (Contribution(b.I, g * b.V),
                Contribution(b.I, white_noise(4 * kb * temp * g)))  # noqa


class _SwitchFlickerHdl(_SwitchHdl):
    """`_SwitchHdl` plus a CONSTANT 1/f current noise `kf/f` on the same
    branch -- white and coloured sources in ONE element under DIFFERENT
    modulations, the case a joint square root got wrong."""
    instparams = _SwitchHdl.instparams + [
        _HdlParameter(name='kf', desc='Flicker PSD at 1 Hz', unit='A^2/Hz',
                      default=0.0)]

    @staticmethod
    def analog(p, n, cp, cn):
        import sympy
        from pycircuit.circuit.hdl import flicker_noise
        b = Branch(p, n)
        ctrl = Branch(cp, cn)
        s = (1 + sympy.tanh((ctrl.V - vth) / vs)) / 2          # noqa: F821
        g = goff + (gon - goff) * s                            # noqa: F821
        return (Contribution(b.I, g * b.V),
                Contribution(b.I, white_noise(4 * kb * temp * g)),  # noqa
                Contribution(b.I, flicker_noise(kf, 1)))       # noqa: F821


class _PllMultPd(Behavioural):
    """Multiplier phase detector, `V(out) = k V(a) V(b)` -- the minimal
    feedback that can pin a VCO's phase.  A PFD is deliberately NOT used:
    A6 step 2 asks about the Jacobian's RANK, and a linear PD is the
    smallest thing that closes the loop."""
    terminals = ('ap', 'an', 'bp', 'bn', 'outp', 'outn')
    params_as = 'p'
    instparams = [_HdlParameter(name='k', desc='PD gain', unit='1/V',
                                default=1.0)]

    @staticmethod
    def analog(p, ap, an, bp, bn, outp, outn):
        ba, bb, bo = Branch(ap, an), Branch(bp, bn), Branch(outp, outn)
        return (Contribution(bo.V, p.k * ba.V * bb.V),)


class _NuMult(Behavioural):
    params_as = 'p'
    instparams = [Parameter(name='k', desc='gain', unit='A/V^2', default=1.0)]

    @staticmethod
    def analog(p, outp, outn, a, an, b, bn):
        return Contribution(Branch(outp, outn).I,
                            p.k * Branch(a, an).V * Branch(b, bn).V)


class _MixedSlopeLo(Behavioural):
    """Two 1/f sources of DIFFERENT slope (0.8 and 2.0) on two branches of
    one element, each times V(lo): a component whose exponent differs
    between entries, stating signed amplitudes (the scale factor outside
    the noise call)."""
    params_as = 'p'
    instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]

    @staticmethod
    def analog(p, a, an, b, bn, lo, lon):
        from pycircuit.circuit.hdl import flicker_noise as _fn
        v = Branch(lo, lon).V
        return (Contribution(Branch(a, an).I, v * _fn(p.k, 0.8)),
                Contribution(Branch(b, bn).I, v * _fn(p.k, 2.0)))


class _Flicker08(Behavioural):
    """A stationary 1/f^0.8 current source."""
    params_as = 'p'
    instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]

    @staticmethod
    def analog(p, a, an):
        from pycircuit.circuit.hdl import flicker_noise as _fn
        return Contribution(Branch(a, an).I, _fn(p.k, 0.8))


class _Flicker20(Behavioural):
    """A stationary 1/f^2 current source."""
    params_as = 'p'
    instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]

    @staticmethod
    def analog(p, a, an):
        from pycircuit.circuit.hdl import flicker_noise as _fn
        return Contribution(Branch(a, an).I, _fn(p.k, 2.0))


class _NuModNoise(Behavioural):
    params_as = 'p'
    instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]

    @staticmethod
    def analog(p, outp, outn, b, bn):
        return Contribution(Branch(outp, outn).I,
                            white_noise((p.k * Branch(b, bn).V) ** 2))



class _SgnPsdFlicker(Behavioural):
    """1/f noise, PSD-specified: `(k V)^2 / f` -- the sign of `V` is gone."""
    params_as = 'p'
    instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]

    @staticmethod
    def analog(p, outp, outn, b, bn):
        from pycircuit.circuit.hdl import flicker_noise as _fn
        return Contribution(Branch(outp, outn).I,
                            _fn((p.k * Branch(b, bn).V) ** 2, 1))


class _SgnAmpFlicker(Behavioural):
    """The SAME PSD, amplitude-specified: `(k V) * flicker_noise(1)`."""
    params_as = 'p'
    instparams = [Parameter(name='k', desc='scale', unit='', default=1.0)]

    @staticmethod
    def analog(p, outp, outn, b, bn):
        from pycircuit.circuit.hdl import flicker_noise as _fn
        return Contribution(Branch(outp, outn).I,
                            p.k * Branch(b, bn).V * _fn(1, 1))


class _PllPhaseDiv(Behavioural):
    """Smooth phase-domain /N: `V(out) = sin(2 pi V(in) / n)`, no limiter."""
    terminals = ('inp', 'inn', 'outp', 'outn')
    params_as = 'p'
    instparams = [_HdlParameter(name='n', desc='Divide ratio', unit='',
                                default=1.0)]

    @staticmethod
    def analog(p, inp, inn, outp, outn):
        import sympy as _sp
        bi, bo = Branch(inp, inn), Branch(outp, outn)
        return (Contribution(bo.V, _sp.sin(2 * _sp.pi * bi.V / p.n)),)
