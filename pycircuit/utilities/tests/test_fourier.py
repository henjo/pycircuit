import pylab
import numpy as np
from numpy.testing import assert_array_almost_equal, assert_array_equal

from pycircuit.utilities import fourier_analysis

doplot = False

def test_fourier_analysis():
    N = 100
    freq1 = 1e3
    freq2 = 2e3
    a0 = 0.9
    a1 = 0.5
    a2 = 0.7
    ph1 = np.pi/2
    ph2 = -np.pi/3

    np.random.seed(100)
    t = np.concatenate(([0], 
                        np.sort(np.random.rand(N)) / freq1, 
                        [1 / freq1]))

    x = a0 + a1 * np.sin(2 * np.pi * freq1 * t + ph1) + \
        a2 * np.sin(2 * np.pi * freq2 * t + ph2)
    
    if doplot:
        pylab.plot(t,x,'+')
        pylab.show()

    X, freqs = fourier_analysis(t, x, harmonics = [-3, -2, -1, 0, 1, 2, 3])

    freqs_ref = np.array([-3*freq1, -freq2, -freq1, 0, freq1, freq2, 3*freq1])

    Xref = np.array(
        [0, 0.5*a2*np.exp(1j*(np.pi/2-ph2)), 0.5*a1*np.exp(1j*(np.pi/2-ph1)),
         a0, 
         0.5*a1*np.exp(1j*(-np.pi/2+ph1)),0.5*a2*np.exp(1j*(-np.pi/2+ph2)),0])
    
    assert_array_equal(freqs, freqs_ref)

    assert_array_almost_equal(X,Xref, 1)
