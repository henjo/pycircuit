# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.


import numpy as np
from pycircuit.post import Waveform
from pycircuit.post.functions import *

from .test_waveform import testdata1, check_func, check_nonscalar_function

## Test unary functions

def test_db20():
    for testdata in testdata1:
        check_func(db20, lambda x: 20*log10(abs(x)), (testdata,))

def test_db10():
    for testdata in testdata1:
        check_func(db10, lambda x: 10*log10(abs(x)), (testdata,))

def test_phase():
    for testdata in testdata1:
        check_func(phase, lambda x: np.angle(x, deg=True), (testdata,))

def test_deriv():
    check_nonscalar_function(deriv)
    
    n = 100
    f = np.array([1.0, 0.5])
    t = np.linspace(0, 1.0, num = n)

    x = Waveform([f, t], np.array([np.sin(2*np.pi * freq * t) for freq in f]))
    xdot_ref = Waveform([f, t], 
                        np.array([2*np.pi*freq * np.cos(2*np.pi * freq * t)
                                    for freq in f]))

    xdot = deriv(x)
    
    error = abs(xdot - xdot_ref[:,:-1])

    maxerror = Waveform(f, 2 * np.pi * f * 5e-2)

    assert ymax(error) < maxerror


def test_cross():
    """Crossings (scipy's interp1d and brenth, imported in `cross` where they
    are used): a sine's rising and falling zero crossings, the second of
    either kind, each row's of a 2-d waveform, and none."""
    phi = np.arange(0, 4 * np.pi, np.pi / 10) - np.pi / 4
    y = Waveform(phi, np.sin(phi))
    assert abs(cross(y)) < 1e-9
    assert abs(cross(y, crosstype=falling) - np.pi) < 1e-9
    assert abs(cross(y, n=1) - np.pi) < 1e-9
    x1 = [np.pi / 4, np.pi / 2]
    rows = cross(Waveform([x1, phi], np.sin(np.vstack([phi for _ in x1]))))
    assert np.allclose(rows.get_x()[0], x1) and np.allclose(rows.get_y(), 0.0, atol=1e-9)
    assert np.isnan(cross(Waveform([[0, 1, 2, 3]], np.array([1, 1, 1, 1]))))
