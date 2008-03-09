import numpy
from numpy import array,concatenate,alltrue,max,min,log10
import scipy
import scipy.interpolate as interpolate
import types
import operator
import pylab

class ResultSet(object):
    """The ResultCollection class handles a set of results"""

    def getResultNames(self):
        pass
    def getResult(self, name):
        pass

class Result(object):
    """The result class manages results from a circuit simulation.

       A result contains the values of one or more signals. Result can
       handle multi dimensional sweeps but not several sweeps.
    """
    class SignalAccesser(object):
        def __init__(self, result):
            self.__result = result
            self.update()
        def update(self):
            self.__dict__.update(dict((k,None) for k in self.__result.getSignalNames()))
        def __getattribute__(self, attr):
            if not attr in ('__dict__', '_SignalAccesser__result', 'update'):
                return self.__result.getSignal(attr)
            return object.__getattribute__(self, attr)

    def __init__(self):
        self.o = self.SignalAccesser(self)

    def getSweepDimensions(self):
        """Return the number of nested sweeps"""
        pass
    def getSweepValues(self, dimension):
        """Get a numpy array of sweep values from sweep dimension dimension"""
        pass

    def getSignalNames(self):
        """Returns a tuple of available signals in the result"""
        pass
    def getSignal(self, name):
        """Returns a Waveform object if the parameter is swept and a scalar otherwise.
           If name is None it will return a dictionary of all signals keyed by
           the signal names
        """

class Waveform(object):
    """The Waveform class handles swept signals. The sweep can be multi dimensional"""

    def __init__(self, x=array([]), y=array([]),
                 xlabels=None, ylabel=None,
                 xunits=None, yunit=None
                 ):
        
        if type(x) != types.ListType:
            x = [x]
            
        dim = len(x)

        assert(dim == len(y.shape))
               
        self._xlist = x
        self._y = y
        self._dim = dim
        self.xlabels = xlabels
        self.ylabel = ylabel
        self.xunits = xunits
        self.yunit = yunit
        
    def getSweepDimensions(self):
        """Return the number of nested sweeps

        >>> w=Waveform([array([1,2]), array([1,2])], array([[3,4],[4,5]]))
        >>> w.getSweepDimensions()
        2
    
        """
        return self._dim

    def getX(self, dim=None):
        """Get X vector of the given sweep dimension

        If no dimension is given, a list of the sweeps in all dimensions is returned

        >>> w=Waveform([array([1,2]), array([1.5,2.5])], array([[3,4],[4,5]]))
        >>> w.getX(0)
        array([1, 2])
        >>> w.getX(1)
        array([ 1.5,  2.5])
        """

        if dim == None:
            return self._xlist
        else:
            return self._xlist[dim]
    
    def getY(self):
        """Get Y vector or n-dimensional array if sweep dimension > 1"""
        return self._y

    def setX(self,value, dim=0):
        "Set X vector"
        self._xlist[dim] = value
    def setY(self,value):
        "Set Y multi-dimensional array"
        self._y = value

    # Operations on Waveform objects
    def __add__(self, a):
        """Add operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w2=Waveform(array([1,2,3]),array([1,1.5,2]))
        >>> w1+w2
        Waveform([1 2 3],[ 4.   6.5  8. ])
        >>> w1+2.0
        Waveform([1 2 3],[ 5.  7.  8.])

        """
        if isinstance(a, Waveform):
            assert(alltrue(concatenate(self._xlist) == concatenate(a.getX())), "X values must be the same")
            return Waveform(self._xlist, self._y + a._y)
        else:
            return Waveform(self._xlist, self._y + a)
    def __radd__(self, a):
        """Reverse add operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> 2.0+w1
        Waveform([1 2 3],[ 5.  7.  8.])

        """
        return self.__add__(a)
    def __sub__(self, a):
        """Subtract operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w2=Waveform(array([1,2,3]),array([1,1.5,2]))
        >>> w1-w2
        Waveform([1 2 3],[ 2.   3.5  4. ])
        >>> w1-2.0
        Waveform([1 2 3],[ 1.  3.  4.])

        """
        if isinstance(a, Waveform):
            assert(alltrue(concatenate(self._xlist) == \
                                 concatenate(a.getX())), "X values must be the same")
            return Waveform(self._xlist, self._y - a._y)
        else:
            return Waveform(self._xlist, self._y - a)
    def __rsub__(self, a):
        """Reverse subtract operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> 2.0-w1
        Waveform([1 2 3],[-1. -3. -4.])

        """
        if isinstance(a, Waveform):
            assert(reduce(operator.__and__, map(lambda x,y: alltrue(x==y), self._xlist, a._xlist)))
            return Waveform(self._xlist, a._y-self._y)
        else:
            return Waveform(self._xlist, a-self._y)
    def __mul__(self, a):
        """Multiplication operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w2=Waveform(array([1,2,3]),array([1,1.5,2]))
        >>> w1*w2
        Waveform([1 2 3],[  3.    7.5  12. ])
        >>> w1*2.0
        Waveform([1 2 3],[  6.  10.  12.])

        """
        if isinstance(a, Waveform):
            assert(reduce(operator.__and__, map(lambda x,y: alltrue(x==y), self._xlist, a._xlist)))
            return Waveform(self._xlist, self._y * a._y)
        else:
            return Waveform(self._xlist, self._y * a)
    def __rmul__(self, a):
        """Reverse multiplication operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> 2.0*w1
        Waveform([1 2 3],[  6.  10.  12.])

        """
        return self.__mul__(a)

    def __div__(self, a):
        """Division operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w2=Waveform(array([1,2,3]),array([1,1.5,2]))
        >>> w1/w2
        Waveform([1 2 3],[ 3.          3.33333333  3.        ])
        >>> w1/2.0
        Waveform([1 2 3],[ 1.5  2.5  3. ])

        """
        if isinstance(a, Waveform):
            assert(reduce(operator.__and__, map(lambda x,y: alltrue(x==y), self._xlist, a._xlist)))
            return Waveform(self._xlist, self._y/a._y)
        else:
            return Waveform(self._xlist, self._y/a)
    def __rdiv__(self, a):
        """Reverse division operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> 2.0/w1
        Waveform([1 2 3],[ 0.66666667  0.4         0.33333333])

        """
        if isinstance(a, Waveform):
            assert(reduce(operator.__and__, map(lambda x,y: alltrue(x==y), self._xlist, a._xlist)))
            return Waveform(self._xlist, a._y/self._y)
        else:
            return Waveform(self._xlist, a/self._y)

    def __abs__(self):
        """Absolute value operator

        >>> w1=Waveform(array([1,2,3]),array([complex(-1,0),complex(0,1),complex(1,-1)]))
        >>> abs(w1)
        Waveform([1 2 3],[ 1.          1.          1.41421356])

        """
        return Waveform(self._xlist, abs(self._y))

    def ymax(self):
        """Returns the maximum y-value

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w1.ymax()
        6
        """
        return max(self._y)

    def ymin(self):
        """Returns the minimum y-value

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w1.ymin()
        3
        """
        return min(self._y)

    def value(self, x):
        """Returns the interpolated y values at x

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w1.value(1.5)
        array([ 4.])
        """
        if self._dim == 1:
            f=scipy.interpolate.interp1d(self._xlist[0], self._y)
        return f(x)

    # Mathematical functions
    def log10(self):
        """Return 20*log10(x)

        >>> w1=Waveform(array([1,2,3]),array([1.0, 10.0]))
        >>> w1.log10()
        Waveform([1 2 3],[ 0.  1.])

        """
        return Waveform(self._xlist, log10(self._y))
        
    def db20(self):
        """Return 20*log10(x)

        >>> w1=Waveform(array([1,2,3]),array([complex(-1,0),complex(0,1),complex(1,-1)]))
        >>> w1.db20()
        Waveform([1 2 3],[ 0.          0.          3.01029996])

        """
        return 20.0*abs(self).log10()

    # Plotting (wrapper around matplotlib)
    def _plot(self, plotfunc, *args, **kvargs):
        import pylab
        if self.getSweepDimensions() == 1:
            p=plotfunc(self.getX(0), self.getY())
        pylab.xlabel('x')
        pylab.ylabel('y')

    def plot(self, *args, **kvargs): self._plot(pylab.plot, *args, **kvargs)
    def semilogx(self, *args, **kvargs): self._plot(pylab.semilogx, *args, **kvargs)
    def semilogy(self, *args, **kvargs): self._plot(pylab.semilogy, *args, **kvargs)
    def loglog(self, *args, **kvargs): self._plot(pylab.loglog, *args, **kvargs)
    
    def __repr__(self):
        if self._dim > 1:
            xlist = self._xlist
        else:
            xlist = self._xlist[0]
        return self.__class__.__name__ + "(" + str(xlist) + "," + str(self.getY()) + ")"

## Waveform functions
def db20(x):
    return 20.0*log10(abs(x))

def db10(x):
    return 10.0*log10(abs(x))

if __name__ == "__main__":
    import doctest
    doctest.testmod()
