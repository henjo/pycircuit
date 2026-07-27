# -*- coding: utf-8 -*-
# Copyright (c) 2026 Pycircuit Development Team
# See LICENSE for details.
"""Every element's ``G`` must be the Jacobian of its ``i``.

An element may state its conductance matrix twice -- as a hand-written ``G()``
stamp, and implicitly as ``i()``, which a differentiating toolkit turns into the
same matrix by autodiff.  Nothing in the type system makes those agree, and when
they disagree Newton is handed a wrong Jacobian: it may still converge, more
slowly or to a slightly different point, so the error hides.

That is exactly what happened to ``VCVS_limited``.  Its ``i()`` computed
``v_out - f'(v_in)*f(v_in)`` and dropped the gain ``g`` entirely, while ``G()``
stamped ``g*f'(u)`` -- and separately ``func.Tanh.fprime`` returned a hard zero,
so the stamp had no input coupling at all.  Both survived for years because the
test suite only ever checked that simulations converged and produced plausible
waveforms, which they did.

This test compares the two directly, by finite difference.  It is the check that
would have caught it.
"""

import numpy as np
import pytest

from pycircuit.circuit import circuit as circuit_module
from pycircuit.circuit.elements import (SubCircuit, VCVS, VCVS_limited, VCCS,
                                        Diode, R, G as Gelem)
from pycircuit.circuit.toolkit import numeric


def _jacobian_of_i(element, x, eps=1e-7):
    """Finite-difference d i / d x, column by column."""
    n = len(x)
    columns = []
    for k in range(n):
        step = np.zeros(n)
        step[k] = eps
        ip = np.asarray(element.i(x + step), dtype=float)
        im = np.asarray(element.i(x - step), dtype=float)
        columns.append((ip - im) / (2 * eps))
    return np.array(columns).T


def _build(element_factory):
    """Instantiate inside a SubCircuit so instance parameters are resolved."""
    saved = circuit_module.default_toolkit
    circuit_module.default_toolkit = numeric
    try:
        cir = SubCircuit(toolkit=numeric)
        cir['E'] = element_factory()
        cir.update_iparv()
        return cir['E']
    finally:
        circuit_module.default_toolkit = saved


## (name, factory, operating points).  Points deliberately straddle the knee of
## any limiting function, since that is where a wrong derivative shows up.
CASES = [
    ('VCVS_limited',
     lambda: VCVS_limited('inp', 'inn', 'outp', 'outn',
                          g=2.0, level=0.5, offset=0.0),
     [[0.1, 0.0, 0.3, 0.0, 0.05],
      [-0.2, 0.1, 0.0, 0.0, 0.0],
      [0.4, 0.0, 1.0, 0.0, 0.1],
      [1.5, 0.0, 0.2, 0.0, 0.0]]),
    ('Diode',
     lambda: Diode('plus', 'minus'),
     [[0.3, 0.0], [0.5, 0.1], [0.0, 0.0]]),
]


@pytest.mark.parametrize('name, factory, points', CASES,
                         ids=[c[0] for c in CASES])
def test_G_is_the_jacobian_of_i(name, factory, points):
    element = _build(factory)
    for point in points:
        x = np.array(point, dtype=float)
        stamped = np.asarray(element.G(x), dtype=float)
        differentiated = _jacobian_of_i(element, x)
        assert np.allclose(stamped, differentiated, atol=1e-5), (
            '%s: G() disagrees with d i/d x at x=%s\n  G()      =\n%s\n'
            '  d i/d x  =\n%s' % (name, point, stamped, differentiated))


def test_vcvs_limited_uses_its_gain():
    """The gain must actually appear in the branch equation.

    Pinned separately because the historical bug was precisely that ``i()``
    ignored ``g``: a test that only compared G against di/dx would have passed
    had *both* dropped it.
    """
    weak = _build(lambda: VCVS_limited('inp', 'inn', 'outp', 'outn',
                                       g=1.0, level=10.0, offset=0.0))
    strong = _build(lambda: VCVS_limited('inp', 'inn', 'outp', 'outn',
                                         g=4.0, level=10.0, offset=0.0))
    x = np.array([0.1, 0.0, 0.0, 0.0, 0.0])
    ## Branch residual is the last entry; with the output held at zero it is
    ## -g*f(v_inn - v_inp), so it must scale with g.
    rw = float(np.asarray(weak.i(x), dtype=float)[-1])
    rs = float(np.asarray(strong.i(x), dtype=float)[-1])
    assert abs(rw) > 1e-12, 'branch residual is identically zero'
    assert abs(rs / rw - 4.0) < 1e-6, (
        'residual did not scale with the gain: g=1 -> %g, g=4 -> %g' % (rw, rs))


def test_tanh_fprime_is_the_derivative_of_f():
    """``func.Tanh.fprime`` must differentiate ``func.Tanh.f``.

    It returned a constant 0 from 2009 until 2026, with the real expression
    unreachable beneath it.
    """
    from pycircuit.circuit import func
    fn = func.Tanh(0.0, 0.5, toolkit=numeric)
    eps = 1e-7
    for x in (-1.0, -0.1, 0.0, 0.1, 1.0):
        fd = (fn.f(x + eps) - fn.f(x - eps)) / (2 * eps)
        assert abs(fn.fprime(x) - fd) < 1e-5, (
            'fprime(%g)=%g but d f/d x=%g' % (x, fn.fprime(x), fd))
