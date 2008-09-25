import numpy as N
from numpy import array,concatenate,alltrue,max,min,log10,arange,pi,sin, \
    sign, where, newaxis, r_, vstack, apply_along_axis, nan, isscalar, rank, inf, sqrt
import scipy
import scipy.interpolate as interpolate
import scipy.optimize as optimize
import types
import operator
import pylab
from copy import copy

# Cartesian operator of  list
def cartesian(listList):
    if listList:
        result = []
        prod = cartesian(listList[:-1])
        for x in prod:
            for y in listList[-1]:
                result.append(x + (y,))
        return result
    return [()]

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

        if dim != len(y.shape):
            raise ValueError("Dimension of x (%s) does not match dimension of y (%s)"%(map(len, x), y.shape))
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

    @property
    def ndim(self):
        """Return the number of nested sweeps"""
        return self._y.ndim

    def getX(self, axis=None):
        """Get X vector of the given sweep dimension

        If no dimension is given, a list of the sweeps in all dimensions is returned

        >>> w=Waveform([array([1,2]), array([1.5,2.5])], array([[3,4],[4,5]]))
        >>> w.getX(0)
        array([1, 2])
        >>> w.getX(1)
        array([ 1.5,  2.5])
        """

        axis = self.getaxis(axis)

        if axis == None:
            return self._xlist
        else:
            return self._xlist[axis]
    
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
    def binaryop(self, op, a):
        """Apply binary operator between self and a"""
        if isinstance(a, Waveform):
            assert(reduce(operator.__and__, map(lambda x,y: alltrue(x==y), self._xlist, a._xlist)))
            return Waveform(self._xlist, op(self._y, a._y), xlabels = self.xlabels)
        else:
            return Waveform(self._xlist, op(self._y, a), xlabels = self.xlabels)

    def rbinaryop(self, op, a):
        """Apply binary operator between a and self"""
        if isinstance(a, Waveform):
            assert(reduce(operator.__and__, map(lambda x,y: alltrue(x==y), self._xlist, a._xlist)))
            return Waveform(self._xlist, op(a._y, self._y), xlabels = self.xlabels)
        else:
            return Waveform(self._xlist, op(a, self._y), xlabels = self.xlabels)

    def __add__(self, a):
        """Add operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w2=Waveform(array([1,2,3]),array([1,1.5,2]))
        >>> w1+w2
        Waveform([1 2 3],[ 4.   6.5  8. ])
        >>> w1+2.0
        Waveform([1 2 3],[ 5.  7.  8.])

        """
        return self.binaryop(operator.__add__, a)
    
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
            ay = a._y
            y = self._y
            if self._y.ndim < a._y.ndim:
                y = self._y[..., newaxis]
            if a._y.ndim < self._y.ndim:
                ay = a._y[..., newaxis]

            return Waveform(self._xlist, y - ay, xlabels = self.xlabels)
        else:
            return Waveform(self._xlist, self._y - a, xlabels = self.xlabels)
    def __rsub__(self, a):
        """Reverse subtract operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> 2.0-w1
        Waveform([1 2 3],[-1. -3. -4.])

        """
        if isinstance(a, Waveform):
            assert(reduce(operator.__and__, map(lambda x,y: alltrue(x==y), self._xlist, a._xlist)))
            return Waveform(self._xlist, a._y-self._y, xlabels = self.xlabels)
        else:
            return Waveform(self._xlist, a-self._y, xlabels = self.xlabels)
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
            return Waveform(self._xlist, self._y * a._y, xlabels = self.xlabels)
        else:
            return Waveform(self._xlist, self._y * a, xlabels = self.xlabels)
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
        return self.binaryop(operator.__div__, a)

    def __rdiv__(self, a):
        """Reverse division operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> 2.0/w1
        Waveform([1 2 3],[ 0.66666667  0.4         0.33333333])

        """
        if isinstance(a, Waveform):
            assert(reduce(operator.__and__, map(lambda x,y: alltrue(x==y), self._xlist, a._xlist)))
            return Waveform(self._xlist, a._y/self._y, xlabels = self.xlabels)
        else:
            return Waveform(self._xlist, a/self._y, xlabels = self.xlabels)

    def __abs__(self):
        """Absolute value operator

        >>> w1=Waveform(array([1,2,3]),array([complex(-1,0),complex(0,1),complex(1,-1)]))
        >>> abs(w1)
        Waveform([1 2 3],[ 1.          1.          1.41421356])

        """
        return Waveform(self._xlist, abs(self._y), xlabels = self.xlabels)

    def __neg__(self):
        """Unary minus operator

        >>> w1=Waveform(array([1,2,3]),array([complex(-1,0),complex(0,1),complex(1,-1)]))
        >>> -w1
        Waveform([1 2 3],[ 1.-0.j -0.-1.j -1.+1.j])

        """
        return Waveform(self._xlist, -self._y, xlabels = self.xlabels)

    def __pow__(self, a):
        """Reverse multiplication operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w1**2
        Waveform([1 2 3],[ 9 25 36])

        """
        return Waveform(self._xlist, self._y ** a, xlabels = self.xlabels)

    def __rpow__(self, a):
        """Reverse multiplication operator

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> 2**w1
        Waveform([1 2 3],[ 8 32 64])

        """
        return Waveform(self._xlist, a ** self._y, xlabels = self.xlabels)

    def __eq__(self, x):
        """Equality operator

        >>> w1=Waveform(array([1,2,3]),array([1,5,2]))
        >>> w2=Waveform(array([1,2,3]),array([1,1.5,2]))
        >>> w1 == w2
        Waveform([1 2 3],[ True False  True])
        >>> w1 == 5
        Waveform([1 2 3],[False  True False])
        
        """
        return self.binaryop(operator.__eq__, x)

    def __lt__(self, x):  return self.binaryop(operator.__lt__, x)
    def __gt__(self, x):  return self.binaryop(operator.__gt__, x)
    def __le__(self, x):  return self.binaryop(operator.__le__, x)
    def __ge__(self, x):  return self.binaryop(operator.__ge__, x)

    def ymax(self, axis=-1):
        """Returns the maximum y-value
        
        Examples
        ========

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w1.ymax()
        6

        >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
        >>> w2.ymax()
        Waveform([1, 2],[6 7])

        """
        return reducedim(self, N.max(self._y, axis=self.getaxis(axis)), axis=self.getaxis(axis))

    def ymin(self, axis=-1):
        """Returns the minimum y-value

        Examples
        ========

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w1.ymin()
        3

        >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
        >>> w2.ymin()
        Waveform([1, 2],[3 4])

        """
        return reducedim(self, N.min(self._y, axis=self.getaxis(axis)))

    def value(self, x, axis = -1, ylabel = None):
        """Returns and interpolated at a given point(s)
        
        `x` can be a number or a waveform with one less than the object
        

        Example
        =======

        1-d waveform

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w1.value(1.5)
        4.0

        2-d waveform
        >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
        >>> w2.value(2.5)
        Waveform([1, 2],[ 4.  5.])
        
        `x` is a waveform
        >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
        >>> w2.value(Waveform([[1, 2]], array([2.5, 3.5])))
        Waveform([1, 2],[ 4.   6.5])

        """
        axis = self.getaxis(axis)
        def findvalue(y):
            return scipy.interpolate.interp1d(self._xlist[axis], y)(x)

        def findvalue_mdimindex(y, i):
            xindex = list(i)
            del xindex[axis]
            xindex = tuple(xindex)
            return scipy.interpolate.interp1d(self._xlist[axis], y)(x._y[xindex])

        if iswave(x):
            newyshape = list(self._y.shape)
            del newyshape[axis]
            newyshape = tuple(newyshape)
            newy = apply_along_axis_with_idx(findvalue_mdimindex,
                                             axis,
                                             self._y).reshape(newyshape)
            return reducedim(self, newy, axis=axis)

        if ylabel == None:
            ylabel = self.ylabel

        outw = applyfunc_and_reducedim(findvalue, self, axis = axis)
        if outw and not isscalar(outw):
            outw.ylabel = ylabel
        return outw

    # Mathematical functions
    def real(self):
        outw = copy(self)
        outw._y = N.real(self._y)
        return outw

    def imag(self):
        outw = copy(self)
        outw._y = N.imag(self._y)
        return outw

    def conjugate(self):
        outw = copy(self)
        outw._y = N.conjugate(self._y)
        return outw

    def log10(self):
        """Return 20*log10(x)

        >>> w1=Waveform(array([1,2,3]),array([1.0, 10.0]))
        >>> w1.log10()
        Waveform([1 2 3],[ 0.  1.])

        """
        return Waveform(self._xlist, log10(self._y), xlabels = self.xlabels)

    def sqrt(self):
        """Return 20*log10(x)

        >>> w1=Waveform(array([1,2,3]),array([1.0, 9.0]))
        >>> w1.sqrt()
        Waveform([1 2 3],[ 1.0  3.])

        """
        return Waveform(self._xlist, sqrt(self._y), xlabels = self.xlabels)
        
    def db20(self):
        """Return x in dB where x is assumed to be a non-power quantity

        >>> w1=Waveform(array([1,2,3]),array([complex(-1,0),complex(0,1),complex(1,-1)]), ylabel='x')
        >>> w1.db20()
        Waveform([1 2 3],[ 0.          0.          3.01029996])
        >>> w1.db20().ylabel
        'db20(x)'
        
        """
        result = 20.0*abs(self).log10()
        if self.ylabel:
            result.ylabel = 'db20(%s)'%self.ylabel
        return result

    def clip(self, xfrom, xto = None, axis=-1):
        """Restrict the waveform to the range defined by xfrom and xto

        >>> w1 = Waveform(array([1.,2.,3.]), array([8., 6., 1.]), ylabel='a')
        >>> w1.clip(2,3)
        Waveform([ 2.  3.],[ 6.  1.])
        >>> w1.clip(1.5, 3)
        Waveform([ 1.5, 2.  3.],[ 7., 6.  1.])

        
        """
        ifrom = self._xlist[axis].searchsorted(xfrom)

        if xto:
            ito = self._xlist[axis].searchsorted(xto)
        else:
            ito = -1

        xaddleft = []
        xaddright = []
        if self._xlist[axis][ifrom] != xfrom:
            pass
        if self._xlist[axis][ito] != xto:
            pass

        newxlist = copy(self._xlist)

        newxlist[axis] = newxlist[axis][ifrom:ito+1]

        ## FIXME, this assumes axis=0
        newy = self._y[ifrom:ito+1]

        return Waveform(newxlist, newy, xunits=self.xunits, yunit=self.yunit, xlabels=self.xlabels, ylabel=self.ylabel)

    # Plotting (wrapper around matplotlib)
    def _plot(self, plotfunc, *args, **kvargs):
        import pylab

        pylab.hold(True)
        for i in N.ndindex(*self._y.shape[:-1]):
            label = ','.join([self.xlabels[axis] + '=' + str(self._xlist[axis][ix]) for axis, ix in enumerate(i)])

            # Limit infinite values
            y = self.getY()[i]
            y[where(y == inf)] = 1e20
            y[where(y == -inf)] = -1e20
            
            p=plotfunc(self.getX(-1), y, label = label, **kvargs)
        pylab.hold(False)
        
        pylab.xlabel(self.xlabels[-1])
        pylab.ylabel(self.ylabel)

    def plot(self, *args, **kvargs): self._plot(pylab.plot, *args, **kvargs)
    def semilogx(self, *args, **kvargs): self._plot(pylab.semilogx, *args, **kvargs)
    def semilogy(self, *args, **kvargs): self._plot(pylab.semilogy, *args, **kvargs)
    def loglog(self, *args, **kvargs): self._plot(pylab.loglog, *args, **kvargs)
    
    @property
    def astable(self):
        """Return a table in text format

        >>> print Waveform(array([1,2,3]),array([3,4,5])).astable
        ====  ===
          x0    y
        ====  ===
           1    3
           2    4
           3    5
        ====  ===

        >>> print Waveform(array([1,2]),array([3,4]), xlabels = ('X',), ylabel = 'Y').astable
        ===  ===
          X    Y
        ===  ===
          1    3
          2    4
        ===  ===

        >>> t=Waveform(array([1,2]),array([3,4]), xlabels = ['X'], ylabel = 'Y').astable

        
        """
        return astable(self)

    def getaxis(self, axis):
        """Look up axis index by name of xlabel names"""
        if type(axis) is types.StringType:
            if axis not in self.xlabels:
                raise Exception('No axis with xlabel %s (%s)'%(axis, str(self.xlabels)))
            return list(self.xlabels).index(axis)
        else:
            return axis

    def reducedimension(self, axes):
        """Reduce given axes by selecting the first element

        >>> w = Waveform([[1,2],[3,4]], array([[1,2],[3,4]]))
        >>> w.reducedimension([0])
        Waveform([[3, 4]],[1 2])
        
        """
        axes = [self.getaxis(axis) for axis in axes]
        
        theslice = list(N.index_exp[:] * self.ndim)
        for axis in axes: theslice[axis] = 0

        w = copy(self)
        w._y = self._y[theslice]

        newxlist = []
        newxlabels = []
        newxunits = []
        for axis in range(self.ndim):
            if axis not in axes:
                newxlist.append(w._xlist[axis])
                if w._xunits:
                    newxunits.append(w._xunits[axis])
                if w._xlabels:
                    newxlabels.append(w._xlabels[axis])

        w._xlist = newxlist
        if w._xunits:
            w._xunits = newxunits
        if w._xlabels:
            w._xlabels = newxlabels
                
        return w
    
    def get_xunits(self):
        if self._xunits != None:
            return self._xunits
        else:
            return len(self._xlist) * ('',)
    def set_xunits(self, units):
        self._xunits = self.__checklabels(units)

    def get_yunit(self):
        if self._yunit != None:
            return self._yunit
        else:
            return ''
    def set_yunit(self, s):
        if type(s) is not types.StringType and s != None:
            raise ValueError('Unit must be a string')
        self._yunit = s

    def get_xlabels(self):
        if self._xlabels != None:
            return self._xlabels
        else:
            return ['x%d'%i for i in range(len(self._xlist))]
    def set_xlabels(self, labels):
        self._xlabels = self.__checklabels(labels)

    def get_ylabel(self):
        if self._ylabel != None:
            return self._ylabel
        else:
            return 'y'
    def set_ylabel(self, s):
        if type(s) is not types.StringType and s != None:
            raise ValueError('Label must be a string')
        self._ylabel = s

    x = property(getX, setX, doc = 'x values')
    y = property(getY, setY, doc = 'y values')
    xlabels = property(get_xlabels, set_xlabels, doc = 'x-axis list of labels for each dimension')
    ylabel = property(get_ylabel, set_ylabel, doc = 'y-axis label')
    xunits = property(get_xunits, set_xunits, doc = 'x-axis list of units for each dimension')
    yunit = property(get_yunit, set_yunit, doc = 'y-axis unit')

    def __getitem__(self, key):
        """Return the value of the given key if the data type is dictionary

        >>> wave = Waveform([[1,2]], array([{'a':1, 'b':2}, {'a':3, 'b':4}]))
        >>> wave['a']
        Waveform([1, 2],[1 3])

        """
        def getitem(d):
            return d[key]
        
        ufuncgetitem = N.vectorize(getitem)
        return Waveform(self._xlist, ufuncgetitem(self._y), xlabels = self.xlabels)

    def __repr__(self):
        if self._dim > 1:
            xlist = self._xlist
        else:
            xlist = self._xlist[0]
        return self.__class__.__name__ + "(" + str(xlist) + "," + str(self.getY()) + ")"

    def __checklabels(self, labels):
        if not labels == None:
            try:
                labels = tuple(labels)
            except:
                raise ValueError('Cannot convert labels to tuples')
            if len(labels) != self._dim:
                raise ValueError('Label list should have the same length (%d) as the number of dimensions (%d)'%(len(labels), self._dim))
            for label in labels:
                if type(label) != types.StringType:
                    raise ValueError('Labels should be of type string')
        return labels
                            

## Waveform functions
def iswave(w):
    """Returns true if argument is a waveform"""
    return isinstance(w, Waveform)

def db20(w):
    return w.db20()
def db10(w):
    return w.db10()
def ymax(w, axis=-1):
    """Returns the maximum y-value

    Examples
    ========

    >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
    >>> ymax(w1)
    6

    >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
    >>> ymax(w2)
    Waveform([1, 2],[6 7])
    """
    return w.ymax(axis=axis)

def ymin(w, axis=-1):
    """Returns the minimum y-value along the given axis

    Examples
    ========

    >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
    >>> ymin(w1)
    3

    >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
    >>> ymin(w2)
    Waveform([1, 2],[3 4])

    """
    
    return w.ymin(axis=axis)

def average(w, axis=-1):
    """Calculate average
    
    Example:

    >>> w1=Waveform([range(2), range(2)],array([[1.0, 3.0], [0.0, 5.0]]))
    >>> average(w1)
    Waveform([0, 1],[ 2.   2.5])

    >>> w1=Waveform([range(2), range(2)],array([[1.0, 3.0], [0.0, 5.0]]), xlabels=['row','col'])
    >>> average(w1, axis='row')
    Waveform([0, 1],[ 0.5  4. ])

    """
    return reducedim(w, N.mean(w._y, axis=w.getaxis(axis)), axis=w.getaxis(axis))

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
    return reducedim(w, N.std(w._y, axis=w.getaxis(axis)), axis=w.getaxis(axis))


raising = 1
falling = 2 
either = 3
def cross(w, crossval = 0.0, n=0, crosstype=either):
    """Calculates the x-axis value where a particular crossing with the specified edge type occurs

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
        finterp = scipy.interpolate.interp1d(x, y)
        return optimize.zeros.brenth(finterp, x[iedge], x[iedge+1])
 
    return applyfunc_and_reducedim(findedge, w - crossval, yunit = w.xunits[0], ylabel = w.xlabels[-1])

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
    return applyfunc(N.angle, w) * 180 / pi

def phaseMargin(g):
    """Calculate phase margin of a open-loop gain vs frequency waveform

    >>> w = 2 * pi * N.logspace(3,8,41)
    >>> w1 = -1e6
    >>> H = Waveform(w, 1.5 * (1 / (1 - 1j*w / w1))**2)
    >>> phaseMargin(H)
    110.376360667
    
    """
    f0 = cross(abs(g), 1.0)
    return phase(-g.value(f0))

def bandwidth(w, db = 3.0, type = 'low'):
    """Calculate bandwidth of transfer as function of frequency

    Examples
    ========

    >>> w = 2 * pi * N.logspace(3,8)
    >>> w1 = -1e6
    >>> H = Waveform(w, 1 / (1 - 1j*w / w1))
    >>> bandwidth(H)
    1000896.9666087811
    
    """
    xmin = min(w._xlist[-1])
    w0 = abs(w.value(xmin))
    
    return cross(abs(w), w0*10**(-db/20.0))

def value(w, x):
    if isinstance(w, Waveform):
        return w.value(x)
    else:
        return w
    

def unityGainFrequency(g):
    """Calculate the frequency where the gain is unity
    """
    return cross(abs(g), 1.0)
    
def IIP3(output, input, fund1, fund2, fund0=None):
    """Calculate input referred third order intermodulation intercept point"""
    s = abs(output)
    if fund0 == None:
        gain = value(s/abs(input), fund1)
    else:
        gain = value(s, abs(fund1)) / value(abs(input), abs(abs(fund1)+fund0))
    return sqrt(s.value(abs(fund1)) * value(s,abs(fund2))**2 / s.value(fund1+2*fund2))/gain

def clip(w, xfrom, xto=None):
    if isinstance(w, Waveform):
        return w.clip(xfrom, xto)
    else:
        return w
    
def applyfunc(func, w):
    if iswave(w):
        outw  = copy(w)
        outw._y = func(outw._y)
        if w.ylabel:
            outw.ylabel = func.__name__ + '(' + w.ylabel + ')'
        return outw
    else:
        return func(w)

def applyfunc_and_reducedim(func, w, axis = -1, ylabel = None, yunit = None):
    """Apply a function that reduces the dimension by one and return a new waveform or float if zero-rank"""

    newyshape = list(w._y.shape)
    del newyshape[axis]
    newy = apply_along_axis(func, axis, w._y).reshape(newyshape)

    if ylabel != None:
        ylabel = func.__name__ + '(' + ylabel + ')'
    
    return reducedim(w, newy, axis=axis, ylabel=ylabel, yunit=yunit)

def reducedim(w, newy, axis=-1, ylabel=None, yunit=None):
    """Reduce the dimension by one and return a new waveform or float if zero-rank"""

    if rank(newy) == 0:
        return N.asscalar(newy)

    if ylabel == None:
        ylabel = w.ylabel

    if yunit == None:
        yunit = w.yunit

    newxlist = list(w._xlist)
    del(newxlist[axis])

    newxlabels = list(w.xlabels)
    del(newxlabels[axis])

    newxunits = list(w.xunits)
    del(newxunits[axis])
        
    return Waveform(newxlist, newy, xlabels = newxlabels, ylabel = ylabel, 
                    xunits = newxunits, yunit = yunit)

def astable(*waveforms):
    """Return a table of one or more waveforms with the same sweeps in text format

    Examples
    ========
    
    >>> w1 = Waveform([range(2)], array([3,4]), ylabel='V1')
    >>> w2 = Waveform([range(2)], array([4,6]), ylabel='V2')
    >>> print astable(w1,w2)
    ====  ====  ====
      x0    V1    V2
    ====  ====  ====
       0     3     4
       1     4     6
    ====  ====  ====
    

    """
    import rsttable

    xvalues = cartesian(waveforms[0]._xlist)
    yvalues = zip(*[list(w._y.flat) for w in waveforms])

    xlabels = waveforms[0].xlabels
    ylabels = [w.ylabel for w in waveforms]
    xunits = waveforms[0].xunits
    yunits = [w.yunit for w in waveforms]

    hasunits = not reduce(operator.__and__, [yunit == '' for yunit in yunits])
    
    if hasunits:
        return rsttable.toRSTtable(map(lambda x,y: list(x) + list(y), 
                                       [xlabels] + [xunits] + xvalues, 
                                       [ylabels] + [yunits] + yvalues), headerrows = 2)
    else:
        return rsttable.toRSTtable(map(lambda x,y: list(x) + list(y), 
                                       [xlabels] + xvalues, 
                                       [ylabels] + yvalues))
    

def apply_along_axis_with_idx(func1d,axis,arr,*args):
    """ Execute func1d(arr[i], i, *args) where func1d takes 1-D arrays
        and arr is an N-d array.  i varies so as to apply the function
        along the given axis for each 1-d subarray in arr.
    """
    arr = N.asarray(arr)
    nd = arr.ndim
    if axis < 0:
        axis += nd
    if (axis >= nd):
        raise ValueError("axis must be less than arr.ndim; axis=%d, rank=%d."
            % (axis,nd))
    ind = [0]*(nd-1)
    i = N.zeros(nd,'O')
    indlist = range(nd)
    indlist.remove(axis)
    i[axis] = slice(None,None)
    outshape = N.asarray(arr.shape).take(indlist)
    i.put(indlist, ind)
    res = func1d(arr[tuple(i.tolist())], tuple(i.tolist()), *args)
    #  if res is a number, then we have a smaller output array
    if isscalar(res):
        outarr = N.zeros(outshape,N.asarray(res).dtype)
        outarr[tuple(ind)] = res
        Ntot = N.product(outshape)
        k = 1
        while k < Ntot:
            # increment the index
            ind[-1] += 1
            n = -1
            while (ind[n] >= outshape[n]) and (n > (1-nd)):
                ind[n-1] += 1
                ind[n] = 0
                n -= 1
            i.put(indlist,ind)
            res = func1d(arr[tuple(i.tolist())], tuple(i.tolist()), *args)
            outarr[tuple(ind)] = res
            k += 1
        return outarr
    else:
        Ntot = N.product(outshape)
        holdshape = outshape
        outshape = list(arr.shape)
        outshape[axis] = len(res)
        outarr = N.zeros(outshape,N.asarray(res).dtype)
        outarr[tuple(i.tolist())] = res
        k = 1
        while k < Ntot:
            # increment the index
            ind[-1] += 1
            n = -1
            while (ind[n] >= holdshape[n]) and (n > (1-nd)):
                ind[n-1] += 1
                ind[n] = 0
                n -= 1
            i.put(indlist, ind)
            res = func1d(arr[tuple(i.tolist())], tuple(i.tolist()), *args)
            outarr[tuple(i.tolist())] = res
            k += 1
        return outarr

def compatible(*args):
    """Return True if the given waveforms have the same x-values

    Example:
    >>> compatible(Waveform(array([1,2,3]),array([3,5,6])), Waveform(array([1,2,3]),array([1,1.5,2])))
    True
    >>> compatible(Waveform(array([1,2]),array([3,5,6])), Waveform(array([1,2,3]),array([1,1.5,2])))
    False
    
    """
    return True

def compose(wlist, x = None, xlabel = None):
    """Compose list of waveforms into a new waveform where the waveform list becomes the outer sweep

    Example:
    >>> wlist=[Waveform(array([1,2,3]),array([3,5,6])), Waveform(array([1,2,3]),array([1,1.5,2]))]
    >>> w = compose(wlist, x = array([1,2]), xlabel = 'index')
    >>> w
    Waveform([array([1, 2]), array([1, 2, 3])],[[ 3.   5.   6. ]
     [ 1.   1.5  2. ]])
    >>> w.xlabels
    ('index', 'x0')
    
    """
    if not compatible(*wlist):
        return ValueError('Waveforms in wlist are not compatible')
    if x != None and len(wlist) != len(x):
        return ValueError('Number of x-values must be the same as the number of waveforms')

    newy = N.array([w.y for w in wlist])

    if x == None:
        newx = [range(len(wlist))] + wlist[0].x
    else:
        newx = [x] + wlist[0].x

    if xlabel == None:
        xlabel = 'composite index'

    return Waveform(newx, newy,
                    xlabels = [xlabel] + list(wlist[0].xlabels), ylabel = wlist[0].ylabel,
                    xunits = [''] + list(wlist[0].xunits), yunit = wlist[0].yunit)
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
