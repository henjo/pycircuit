"""This module contains functions that operates on wave objects or scalars"""

from waveform import Waveform, reducedim, applyfunc, applyfunc_and_reducedim,\
    iswave, wavefunc
import numpy as np
from numpy import array, pi, sign, alltrue, where, arange, vstack, \
    sin, log10, sqrt, nan
import scipy as sp
import scipy.optimize as optimize

def db10(w):
    """Return x in dB where x is assumed to be a non-power quantity

    >>> w1=Waveform(array([1,2,3]),array([complex(-1,0),complex(0,1),2]), ylabel='x')
    >>> db10(w1)
    Waveform([1 2 3],[ 0.          0.          3.01029996])
    >>> db10(w1).ylabel
    'db10(x)'
        
    """
    return applyfunc(lambda x: 10.0*log10(abs(x)), w, 'db10')

def db20(w):
    """Return x in dB where x is assumed to be a non-power quantity

    >>> w1=Waveform(array([1,2,3]),array([complex(-1,0),complex(0,1),complex(1,-1)]), ylabel='x')
    >>> db20(w1)
    Waveform([1 2 3],[ 0.          0.          3.01029996])
    >>> db20(w1).ylabel
    'db20(x)'
        
    """
    return applyfunc(lambda x: 20.0*log10(abs(x)), w, 'db20')

@wavefunc
def ymax(w, axis=-1):
    return w

@wavefunc
def ymin(w, axis=-1):
    return w

@wavefunc
def value(w, x):
    return w

raising = 1
falling = 2
either = 3
def cross(w, crossval = 0.0, n=0, crosstype=either):
    """Calculates the x-axis value where a particular crossing with the
    specified edge type occurs

    Examples:

    1-d waveform

    >>> phi = arange(0, 4*pi, pi/10)-pi/4
    >>> y = Waveform(phi, sin(phi))
    >>> cross(y)
    0.0
    >>> cross(y, crosstype=falling)
    3.1415926535897931

    2-d waveform

    >>> x1 = [pi/4,pi/2]
    >>> x2 = arange(0, 4*pi, pi/10)-pi/4
    >>> phi = vstack([x2 for p in x1])
    >>> y = Waveform([x1,x2], sin(phi))
    >>> cross(y)
    Waveform([0.78539816339744828, 1.5707963267948966],[ 0.  0.])

    No crossing

    >>> cross(Waveform([[0,1,2,3]], array([1,1,1,1])))
    nan

    TODO: handle case where x-values are exactly at the crossing

    """

    x = w._xlist[-1]

    def findedge(y):
        ## Find edges
        if crosstype == either:
            edges = sign(y[:-1]) != sign(y[1:])
        elif crosstype == raising:
            edges = sign(y[:-1]) < sign(y[1:])
        elif crosstype == falling:
            edges = sign(y[:-1]) > sign(y[1:])

        if alltrue(edges == False):
            return nan

        iedge = where(edges)[0][n]

        ## Find exact location of the crossing using interpolated x-values
        finterp = sp.interpolate.interp1d(x, y)
        return optimize.zeros.brenth(finterp, x[iedge], x[iedge+1])

    return applyfunc_and_reducedim(findedge, w - crossval, yunit = w.xunits[0],
                                   ylabel = w.xlabels[-1])

def phase(w):
    """Return argument in degrees of complex values

    Example:

    >>> phase(1)
    0.0
    >>> phase(complex(0,1))
    90.0
    >>> phase(Waveform((range(3),), array([1, complex(1,1), complex(0,-1)])))
    Waveform(([0, 1, 2],),[  0.  45. -90.])

    """
    return applyfunc(np.angle, w) * 180 / pi

def phaseMargin(g):
    """Calculate phase margin of a open-loop gain vs frequency waveform

    >>> w = 2 * pi * np.logspace(3,8,41)
    >>> w1 = -1e6
    >>> H = Waveform(w, 1.5 * (1 / (1 - 1j*w / w1))**2)
    >>> '%0.4g'%phaseMargin(H)
    '110.4'

    """
    f0 = cross(abs(g), 1.0)
    return phase(-g.value(f0))

def bandwidth(w, db = 3.0, type = 'low'):
    """Calculate bandwidth of transfer as function of frequency

    Examples
    ========

    >>> w = 2 * pi * np.logspace(3,8)
    >>> w1 = -1e6
    >>> H = Waveform(w, 1 / (1 - 1j*w / w1))
    >>> bandwidth(H)
    1000896.9666087811

    """
    xmin = min(w._xlist[-1])
    w0 = abs(w.value(xmin))

    return cross(abs(w), w0*10**(-db/20.0))


def unityGainFrequency(g):
    """Calculate the frequency where the gain is unity
    """
    return cross(abs(g), 1.0)

def IM3(w, fund1, fund2, fund0=None):
    """Return input referred third order intermodulation tone

    The intermodulation product is evaluated at fund1 + 2 * fund2

    """
    return value(abs(w), fund1 + 2 * fund2)

def IM2(w, fund1, fund2, fund0=None):
    """Return input referred third order intermodulation tone

    The intermodulation product is evaluated at fund1 + fund2

    """
    return value(abs(w), fund1 + fund2)

def IIP3(output, input, fund1, fund2, fund0=None):
    """Calculate input referred third order intermodulation intercept point

    The intermodulation product is evaluated at fund1 + 2 * fund2

    """
    s = abs(output)
    if fund0 == None:
        gain = value(s/abs(input), fund1)
    else:
        gain = value(s, abs(fund1)) / value(abs(input), abs(abs(fund1)+fund0))
    return sqrt(s.value(abs(fund1)) * value(s,abs(fund2))**2 /
                value(s, fund1 + 2 * fund2)) / gain

def IIP2(output, input, fund1, fund2, fund0=None):
    """Calculate input referred second order intermodulation intercept point

    The intermodulation product is evaluated at fund1 + fund2

    """
    s = abs(output)
    if fund0 == None:
        gain = value(s/abs(input), fund1)
    else:
        gain = value(s, abs(fund1)) / value(abs(input), abs(abs(fund1)+fund0))

    return value(s, abs(fund1)) * value(s, abs(fund2)) \
           / value(s, fund1 + fund2) / gain

def clip(w, xfrom, xto=None):
    if isinstance(w, Waveform):
        return w.clip(xfrom, xto)
    else:
        return w

def average(w, axis=-1):
    """Calculate average 

    Example:

    >>> w1=Waveform([range(2), range(2)],array([[1.0, 3.0], [0.0, 5.0]]))
    >>> average(w1)
    Waveform([0, 1],[ 2.   2.5])

    >>> w1=Waveform([range(2), range(2)],array([[1.0, 3.0], [0.0, 5.0]]), \
                    xlabels=['row','col'])
    >>> average(w1, axis='row')
    Waveform([0, 1],[ 0.5  4. ])

    """
    return reducedim(w, np.mean(w._y, axis=w.getaxis(axis)),
                     axis=w.getaxis(axis))

def stddev(w, axis=-1):
    """Calculate the standard deviation

    Returns the standard deviation over the highest dimension, a measure of the
    spread of a distribution.

    Examples
    ========

    >>> w1=Waveform([range(2), range(4)], array([[1,2,3,4],[1,1,1,1]]))
    >>> stddev(w1)
    Waveform([0, 1],[ 1.11803399  0.        ])


    """
    return reducedim(w, np.std(w._y, axis=w.getaxis(axis)),
                     axis=w.getaxis(axis))


def deriv(w):
    """Calculate derivative of a waveform with respect to the inner x-axis"""
    return w.deriv()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
