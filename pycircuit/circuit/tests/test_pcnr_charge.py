# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""PCNR with a charge-storing participant (hdl.md sec. 9, phase A1).

A junction device with capacitance is the normal case, not an exotic
one, and until now `pcnr.augmented_system` refused it outright.  It is
allowed because the resulting Newton system is exactly consistent:

* `g_MNA`'s charge part is a function of ``x_MNA`` and its derivative
  ``Geq`` stays in ``J_MNA/MNA``;
* the device's resistive current is a function of ``v_lim`` alone and
  its derivative goes to ``J_MNA/lim``;
* ``J_lim/lim`` is the identity, ``J_lim/MNA`` the incidence row.

This module proves that by finite-differencing the whole augmented
residual, and then checks the thing a user cares about: PCNR on and off
produce the same transient.
"""

import warnings

import numpy as np
import pytest
import sympy
from numpy.testing import assert_allclose

import pycircuit.circuit.circuit
from pycircuit.circuit.toolkit import numeric
from pycircuit.circuit import gnd
from pycircuit.circuit.dcanalysis import DC
from pycircuit.circuit.transient import Transient
from pycircuit.circuit.elements import SubCircuit, VS, VSin, R
from pycircuit.circuit import pcnr as _pcnr
from pycircuit.circuit.hdl import (Behavioural, Branch, Contribution, ddt,
                                   vt)
from pycircuit.utilities.param import Parameter


class DiodeWithCj(Behavioural):
    """A diode that also stores charge -- the case A1 unblocks."""
    instparams = [Parameter(name='IS', desc='Sat current', unit='A',
                            default=1e-13),
                  Parameter(name='cj', desc='Junction capacitance', unit='F',
                            default=1e-12)]

    @staticmethod
    def analog(plus, minus):
        b = Branch(plus, minus)
        return Contribution(b.I,
                            IS * (sympy.exp(b.V / vt()) - 1)   # noqa: F821
                            + ddt(cj * b.V))                   # noqa: F821


def _circuit(vsrc=1.0, sinusoid=False):
    pycircuit.circuit.circuit.default_toolkit = numeric
    c = SubCircuit()
    na, nb = c.add_node('a'), c.add_node('b')
    if sinusoid:
        c['vs'] = VSin(na, gnd, va=vsrc, freq=1e6)
    else:
        c['vs'] = VS(na, gnd, v=vsrc)
    c['R'] = R(na, nb, r=1e3)
    c['D'] = DiodeWithCj(nb, gnd, IS=1e-13, cj=1e-11)
    c.update_iparv()
    return c


def test_charge_storing_element_participates():
    """It qualifies now -- and it really does store charge, so this is
    the case the layer used to refuse."""
    assert DiodeWithCj.pcnr_junctions == ((0, 1),)
    c = _circuit()
    js = _pcnr.pcnr_junctions(c)
    assert len(js) == 1 and js[0][0] == 'D'
    Cm = np.asarray(c['D'].C(np.zeros(2)), float)
    assert np.any(np.abs(Cm) > 0.0), Cm


def _augmented(c, x, v_lim, h, q_last):
    """Assemble (g, J) of the augmented system the way the transient
    does: backward Euler companion terms recomputed at THIS x."""
    from pycircuit.circuit.circuit import defaultepar
    epar = defaultepar
    q = np.asarray(c.q(x, epar), float)
    Cm = np.asarray(c.C(x, epar), float)
    iq = (q - q_last) / h
    Geq = Cm / h
    u = np.asarray(c.u(0.0, epar, analysis='tran'), float)
    js = _pcnr.pcnr_junctions(c)
    g_mna, g_lim, J_mm, J_ml, J_lm, _didv = _pcnr.augmented_system(
        c, x, v_lim, js, epar, u_extra=iq + u, dense_blocks=True,
        J_extra=Geq)
    n, k = len(x), len(js)
    g = np.concatenate([np.asarray(g_mna, float), np.asarray(g_lim, float)])
    J = np.zeros((n + k, n + k))
    J[:n, :n] = J_mm
    J[:n, n:] = J_ml
    J[n:, :n] = J_lm
    J[n:, n:] = np.eye(k)
    return g, J


def test_augmented_jacobian_is_the_derivative_with_charge():
    """The claim in pcnr.augmented_system's comment, finite-differenced:
    with a charge-storing participant, J == dg/dx over BOTH blocks --
    the MNA unknowns and the limited unknown."""
    c = _circuit()
    n = c.n
    x = np.array([1.0, 0.6, -4e-4][:n]) if n == 3 else np.linspace(0.2, 0.6, n)
    v_lim = np.array([0.55])
    h = 1e-9
    q_last = np.zeros(n)

    g0, J = _augmented(c, x, v_lim, h, q_last)
    k = len(v_lim)
    eps = 1e-7
    Jfd = np.zeros_like(J)
    for j in range(n):                      # d/dx_MNA
        dx = np.zeros(n); dx[j] = eps
        gp, _ = _augmented(c, x + dx, v_lim, h, q_last)
        gm, _ = _augmented(c, x - dx, v_lim, h, q_last)
        Jfd[:, j] = (gp - gm) / (2 * eps)
    for j in range(k):                      # d/dv_lim
        dv = np.zeros(k); dv[j] = eps
        gp, _ = _augmented(c, x, v_lim + dv, h, q_last)
        gm, _ = _augmented(c, x, v_lim - dv, h, q_last)
        Jfd[:, n + j] = (gp - gm) / (2 * eps)

    scale = max(1.0, float(np.max(np.abs(J))))
    assert_allclose(J, Jfd, rtol=2e-4, atol=1e-6 * scale)


def test_charge_participant_dc_matches_non_pcnr():
    c = _circuit(vsrc=5.0)
    v_off = float(DC(c, toolkit=numeric, pcnr=False).solve().v('b'))
    v_on = float(DC(c, toolkit=numeric, pcnr=True).solve().v('b'))
    ## PCNR's own convergence tolerance, not a modelling difference: the
    ## two solvers stop at slightly different points on the same curve
    ## (the hand-written Diode shows the same ~2e-6 spread).
    assert_allclose(v_on, v_off, rtol=1e-5)
    assert 0.4 < v_off < 0.8, v_off        # a real forward-biased junction


def test_charge_participant_transient_matches_non_pcnr():
    """The user-visible claim: with capacitance present, a PCNR transient
    reproduces the ordinary one.  This raised NotImplementedError before
    phase A1."""
    def run(pcnr):
        c = _circuit(vsrc=2.0, sinusoid=True)
        tran = Transient(c, toolkit=numeric, pcnr=pcnr)
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            res = tran.solve(tend=2e-6, timestep=2e-8, fixed_timestep=True)
        return np.asarray(res.v('b').y, float)

    y_off = run(False)
    y_on = run(True)
    assert y_on.shape == y_off.shape
    assert_allclose(y_on, y_off, rtol=1e-5, atol=1e-9)
    assert np.ptp(y_off) > 0.1             # the waveform actually swings


def test_charge_participant_transient_used_to_be_refused():
    """Guard the refusal's removal: the layer must no longer raise for a
    charge-storing participant (and the reason it may not is that the
    Jacobian test above passes)."""
    c = _circuit()
    from pycircuit.circuit.circuit import defaultepar
    js = _pcnr.pcnr_junctions(c)
    n = c.n
    x = np.zeros(n)
    Cm = np.asarray(c.C(x, defaultepar), float)
    _pcnr.augmented_system(c, x, np.zeros(len(js)), js, defaultepar,
                           u_extra=np.zeros(n), dense_blocks=True,
                           J_extra=Cm / 1e-9)      # J_extra set == transient


def test_limexp_is_what_makes_a_pcnr_participant_robust():
    """A correction to an earlier claim in hdl.md, kept as a test.

    It said that with PCNR on, ``limexp``'s clamp "is simply never
    reached, because PCNR keeps the iterate in range".  Measurement says
    otherwise, and the mechanism is worth pinning:
    ``pcnr.augmented_system`` assembles ``cir.i(x)`` -- which INCLUDES
    the participant -- and then subtracts that device's own ``i(sub)``
    again, because its current is about to be re-stamped at ``v_lim``.
    The cancellation is exact in arithmetic and worthless in floating
    point once the term is ``inf``: ``inf - inf = nan`` poisons the whole
    system.  PCNR bounds the LIMITED quantity; it does not bound the node
    voltage at which the device's own ``i()`` is evaluated during
    assembly.  ``limexp`` does.
    """
    from pycircuit.circuit.hdl import limexp

    def diode(expfn):
        class _D(Behavioural):
            instparams = [Parameter(name='IS', desc='Is', unit='A',
                                    default=1e-13)]

            @staticmethod
            def analog(plus, minus):
                b = Branch(plus, minus)
                return Contribution(b.I,
                                    IS * (expfn(b.V / vt()) - 1))  # noqa
        return _D

    pycircuit.circuit.circuit.default_toolkit = numeric
    raw, lim = diode(sympy.exp)('p', 'n'), diode(limexp)('p', 'n')
    raw.update_iparv(); lim.update_iparv()
    x20 = np.array([20.0, 0.0])
    assert not np.isfinite(np.asarray(raw.i(x20), float)[0])
    assert np.all(np.isfinite(np.asarray(lim.i(x20), float)))

    def solve(cls, pcnr):
        c = SubCircuit()
        na, nb = c.add_node('a'), c.add_node('b')
        c['vs'] = VS(na, gnd, v=20.0)      # hard forward drive
        c['R'] = R(na, nb, r=1.0)
        c['D'] = cls(nb, gnd, IS=1e-13)
        c.update_iparv()
        return float(DC(c, toolkit=numeric, pcnr=pcnr).solve().v('b'))

    ## Both models solve fine WITHOUT pcnr -- the ordinary path's
    ## continuation ladders handle them.
    v_ref = solve(diode(limexp), False)
    assert 0.7 < v_ref < 1.0, v_ref

    ## With pcnr, the raw-exp model dies in assembly; limexp survives and
    ## lands on the same answer.
    with pytest.raises(Exception):
        solve(diode(sympy.exp), True)
    assert_allclose(solve(diode(limexp), True), v_ref, rtol=1e-5)
