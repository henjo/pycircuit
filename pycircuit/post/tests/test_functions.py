# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from nose.tools import *

import numpy as np
from pycircuit.post import Waveform
from pycircuit.post.functions import *

from test_waveform import testdata1, check_func, check_nonscalar_function

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
