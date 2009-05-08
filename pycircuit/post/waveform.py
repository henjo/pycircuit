# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

"""The waveform module contains classes for handling simulation results in 
the form of X-Y data. The classes can handle results from multi-dimensional 
sweeps.
"""

import numpy as np
from numpy import array,concatenate,alltrue,max,min,log10,arange,pi,sin, \
    sign, where, newaxis, r_, vstack, apply_along_axis, nan, isscalar, rank, \
    inf, isscalar
import scipy as sp
import scipy.interpolate as interpolate
import types
import operator
import pylab
from copy import copy

class Waveform(object):
    """The Waveform class handles swept signals. The sweep can be multi 
    dimensional.

    The waveform class can handle both n-dimensional gridded data or data 
    where the inner dimension has variable length.

    Examples:

    Initiating N-dimensional gridded data:

    >>> w = Waveform([array([1,2]), array([3,4])], array([[3,4],[4,5]]))

    Initiating N-dimensional data where the inner dimension has variable length

    >>> x0 = array([array([1,1]), array([2])], dtype=object)
    >>> x1 = array([array([1,2]), array([1])], dtype=object)
    >>> w = Waveform([x0, x1], array([array([3,4]),array([4])], dtype = object))
    
    """

    def __init__(self, x=array([],), y=array([]),
                 xlabels=None, ylabel=None,
                 xunits=None, yunit=None):
        
        if type(x) is np.ndarray:
            x = [x]
            
        if type(y) is not np.ndarray:
            y = array(y)

        xlist = [array(xelement) for xelement in x]

        self.ragged = (y.dtype == object) and  \
            set([xe.shape for xe in xlist + [y]]) == set([y.shape])
        
        dim = len(x)

        if not self.ragged and dim != len(y.shape):
            raise ValueError("Dimension of x (%s) does not match dimension of"
                             " y (%s)"%(map(len, x), y.shape))
        self._xlist = xlist
        self._y = y
        self._dim = dim
        self.xlabels = xlabels
        self.ylabel = ylabel
        self.xunits = xunits
        self.yunit = yunit

    ## numpy array interface

    __array_priority__ = 100.0

    def __array__(self): return self._y
        
    def __array_wrap__(self, arr, context=None):
        return Waveform(list(self._xlist), arr, 
                        xlabels=self.xlabels, ylabel=self.ylabel,
                        xunits=self.xunits, yunit=self.yunit)
        
    @property
    def ndim(self):
        """Return the number of nested sweeps"""
        return self._y.ndim

    @property
    def shape(self): return tuple((len(x) for x in self.x))

    def get_x(self, axis=None):
        """Get X vector of the given sweep dimension

        If no dimension is given, a list of the sweeps in all dimensions 
        is returned

        >>> w=Waveform([array([1,2]), array([1.5,2.5])], array([[3,4],[4,5]]))
        >>> w.get_x(0)
        array([1, 2])
        >>> w.get_x(1)
        array([ 1.5,  2.5])
        """

        if axis == None:
            return self._xlist
        else:
            axis = self.getaxis(axis)
            return self._xlist[axis]
    
    def get_y(self):
        """Get Y vector or n-dimensional array if sweep dimension > 1"""
        return self._y

    def set_x(self,value, axis=-1):
        "Set X vector"

        ## Swap order if new xvalues are falling
        if value[0] > value [-1]:
            self._xlist[axis] = value[-1::-1]

            if axis != -1:
                raise Exception("Can only swap order if axis=-1")

            self._y = self._y[..., -1::-1]            

    def set_y(self,value):
        "Set Y multi-dimensional array"
        self._y = value

    def map_x(self, func, axis=-1, xlabel = ''):
        """Apply function func on x-values of the given axis"""
        newxlist = copy(self._xlist)
        newxlist[axis] = func(newxlist[axis])

        newxlabels = list(self._xlabels)
        newxlabels[axis] = xlabel

        return Waveform(newxlist, copy(self._y), xlabels = newxlabels, 
                        yunit = self.yunit, ylabel = self.ylabel)


    ## Operations on Waveform objects
    def binaryop(self, op, a, ylabel = None, yunit = None, reverse = False):
        """Apply binary operator between self and a"""
        if isinstance(a, Waveform):
            assert(reduce(operator.__and__, 
                          map(lambda x,y: alltrue(x==y), 
                              self._xlist, a._xlist))), \
                "x-axes of the arguments must be the same"
            ay = a._y
        else:
            ay = a 

        if reverse:
            result = op(ay, self._y)
        else:
            result = op(self._y, ay)            

        return Waveform(self._xlist, result, 
                        xlabels = self.xlabels, xunits = self.xunits,
                        ylabel = ylabel, yunit = yunit)

    ## Unary operators
    def __abs__(self):     
        return Waveform(self._xlist, abs(self._y), xlabels = self.xlabels, 
                        xunits = self.xunits,
                        ylabel = 'abs(%s)'%self.ylabel, yunit = self.yunit)
    def __neg__(self):
        return Waveform(self._xlist, -self._y, xlabels = self.xlabels, 
                        xunits = self.xunits,
                        ylabel = '-%s'%self.ylabel, yunit = self.yunit)
    ## Binary operators
    def __add__(self, a):  return self.binaryop(operator.__add__, a)
    def __radd__(self, a): 
        return self.binaryop(operator.__add__, a, reverse=True)
    def __sub__(self, a):  return self.binaryop(operator.__sub__, a)
    def __rsub__(self, a): 
        return self.binaryop(operator.__sub__, a, reverse=True)
    def __mul__(self, a):  return self.binaryop(operator.__mul__, a)
    def __rmul__(self, a): 
        return self.binaryop(operator.__mul__, a, reverse=True)
    def __div__(self, a):  return self.binaryop(operator.__div__, a)
    def __rdiv__(self, a): 
        return self.binaryop(operator.__div__, a, reverse=True)
    def __pow__(self, a):  return self.binaryop(operator.__pow__, a)
    def __rpow__(self, a): 
        return self.binaryop(operator.__pow__, a, reverse=True)
    def __eq__(self, x):   return self.binaryop(operator.__eq__, x)
    def __lt__(self, x):   return self.binaryop(operator.__lt__, x)
    def __gt__(self, x):   return self.binaryop(operator.__gt__, x)
    def __le__(self, x):   return self.binaryop(operator.__le__, x)
    def __ge__(self, x):   return self.binaryop(operator.__ge__, x)

    def xmax(self, axis=-1):
        """Returns the maximum x-value over one axis
        
        Examples:

        >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
        >>> w2.xmax()
        4

        """
        return np.max(self._xlist[axis])

    def xmin(self, axis=-1):
        """Returns the minimum x-value over one axis
        
        Examples:

        >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
        >>> w2.xmin()
        2

        """
        return np.min(self._xlist[axis])

    def ymax(self, axis=-1):
        """Returns the maximum y-value over one axis
        
        Examples:

        >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
        >>> w2.ymax()
        Waveform(array([1, 2]), array([6, 7]))

        """
        return reducedim(self, np.max(self._y, axis=self.getaxis(axis)), 
                         axis=self.getaxis(axis))

    def ymin(self, axis=-1):
        """Returns the minimum y-value over one axis

        Examples:

        >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
        >>> w2.ymin()
        Waveform(array([1, 2]), array([3, 4]))

        """
        return reducedim(self, np.min(self._y, axis=self.getaxis(axis)), 
                         axis=self.getaxis(axis))

    def value(self, x, axis = -1, ylabel = None):
        """Returns and interpolated at the given x-value
        
        *x* can be a number or a waveform where the number of dimensions of x is
        is one less than the waveform it is operating on
           
        
        Examples:

        1-d waveform

        >>> w1=Waveform(array([1,2,3]),array([3,5,6]))
        >>> w1.value(1.5)
        4.0

        2-d waveform
        >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
        >>> w2.value(2.5)
        Waveform(array([1, 2]), array([ 4.,  5.]))
        
        `x` is a waveform
        >>> w2=Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
        >>> w2.value(Waveform([[1, 2]], array([2.5, 3.5])))
        Waveform(array([1, 2]), array([ 4. ,  6.5]))

        """
        axis = self.getaxis(axis)
        def findvalue(y):
            if len(self._xlist[axis]) == 1 and axis == -1:
                return y[-1]
            res = sp.interpolate.interp1d(self._xlist[axis], y)(x)
            try:
                return np.asscalar(res)
            except TypeError:
                return res

        def findvalue_mdimindex(y, i):
            xindex = list(i)
            del xindex[axis]
            xindex = tuple(xindex)
            res = sp.interpolate.interp1d(self._xlist[axis], y)(x._y[xindex])
            try:
                return np.asscalar(res)
            except TypeError:
                return res
            

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

    def clip(self, xfrom, xto = None, axis=-1):
        """Restrict the waveform to the range defined by xfrom and xto

        >>> w1 = Waveform(array([1.,2.,3.]), array([8., 6., 1.]), ylabel='a')
        >>> w1.clip(2,3)
        Waveform(array([ 2.,  3.]), array([ 6.,  1.]))
        >>> w1.clip(1.5, 3)
        Waveform(array([ 1.5,  2. ,  3. ]), array([ 7.,  6.,  1.]))
        
        """
        ifrom = self._xlist[axis].searchsorted(xfrom)

        if xto:
            ito = self._xlist[axis].searchsorted(xto, side='right')
        else:
            ito = 0

        newxlist = copy(self._xlist)
        newxlist[axis] = newxlist[axis][ifrom:ito]

        newy = self._y[onedim_index(slice(ifrom,ito), axis, self.ndim)]
        
        ## Add new items on left and right side if the clip limit
        ## does not coincidence with an already present x-value
        if self._xlist[axis][ifrom] != xfrom:
            newxlist[axis] = np.insert(newxlist[axis], ifrom-1, xfrom)
            func = lambda y: [np.interp(xfrom, self._xlist[axis], y)]
            yfrom = np.apply_along_axis(func, axis, self._y)
            newy = np.concatenate((yfrom, newy),axis=axis)

        if self._xlist[axis][ito-1] != xto:
            newxlist[axis] = np.insert(newxlist[axis], ito, xto)
            func = lambda y: [np.interp(xto, self._xlist[axis], y)]
            yto = np.apply_along_axis(func, axis, self._y)
            newy = np.concatenate((newy, yto),axis=axis)

        return Waveform(newxlist, newy, xunits=self.xunits, yunit=self.yunit, 
                        xlabels=self.xlabels, ylabel=self.ylabel)

    # Mathematical functions
    def real(self): return np.real(self)
    def imag(self): return np.imag(self)
    def conjugate(self): return np.conjugate(self)

    def deriv(self):
        """Calculate derivative of a waveform with respect to the inner x-axis"""

        newxlist = copy(self._xlist)

        newxlist[-1] = newxlist[-1][:-1]

        dydx = np.diff(self.y, axis=-1) / np.diff(self._xlist[-1])
        
        return Waveform(newxlist, dydx, xlabels = self.xlabels)

    # Plotting (wrapper around matplotlib)
    def _plot(self, plotfunc, *args, **kvargs):
        import pylab

        pylab.hold(True)
        for i in np.ndindex(*self._y.shape[:-1]):
            label = ','.join([self.xlabels[axis] + '=' + \
                    str(self._xlist[axis][ix]) for axis, ix in enumerate(i)])

            # Limit infinite values
            y = self.get_y()[i]
            y[where(y == inf)] = 1e20
            y[where(y == -inf)] = -1e20

#            if 'label' not in kvargs:
#                kvargs['label'] = label
            
            p=plotfunc(self.get_x(-1), y, **kvargs)

        pylab.hold(False)
        
        xlabel = self.xlabels[-1]
        if self.xunits[-1] != '':
            xlabel += ' [%s]'%self.xunits[-1]
        pylab.xlabel(xlabel)
        
        ylabel = self.ylabel
        if self.yunit != '':
            ylabel += ' [%s]'%self.yunit
        pylab.ylabel(ylabel)

    def plot(self, *args, **kvargs): self._plot(pylab.plot, *args, **kvargs)
    def semilogx(self, *args, **kvargs): self._plot(pylab.semilogx, *args, **kvargs)
    def semilogy(self, *args, **kvargs): self._plot(pylab.semilogy, *args, **kvargs)
    def loglog(self, *args, **kvargs): self._plot(pylab.loglog, *args, **kvargs)

    def stem(self, *args, **kvargs): self._plot(pylab.stem, *args, **kvargs)
    
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
        >>> print Waveform(array([1,2]),array([3,4]), xlabels = ('X',), \
                           ylabel = 'Y').astable
        ===  ===
        X    Y  
        ===  ===
          1    3
          2    4
        ===  ===
        >>> t=Waveform(array([1,2]),array([3,4]), xlabels = ['X'], \
                       ylabel = 'Y').astable

        
        """
        return astable(self)

    def getaxis(self, axis):
        """Look up axis index by name of xlabel names"""
        if type(axis) is types.StringType:
            if axis not in self.xlabels:
                raise Exception('No axis with xlabel %s (%s)'%(axis, str(self.xlabels)))
            return list(self.xlabels).index(axis)
        elif type(axis) is types.IntType:
            if axis >= self.ndim:
                raise ValueError('axis %d >= number of dimensions'%axis)
            return axis
        else:
            raise ValueError('Axis %s must be a string or an integer'%str(axis))

    def reducedimension(self, axes):
        """Reduce given axes by selecting the first element

        >>> w = Waveform([[1,2],[3,4]], array([[1,2],[3,4]]))
        >>> w.reducedimension([0])
        Waveform([array([3, 4])], array([1, 2]))
        
        """
        axes = [self.getaxis(axis) for axis in axes]
        
        theslice = list(np.index_exp[:] * self.ndim)
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
            return tuple(['x%d'%i for i in range(len(self._xlist))])
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

    x = property(get_x, set_x, doc = 'x values')
    y = property(get_y, set_y, doc = 'y values')
    xlabels = property(get_xlabels, set_xlabels, \
                           doc = 'x-axis list of labels for each dimension')
    ylabel = property(get_ylabel, set_ylabel, doc = 'y-axis label')
    xunits = property(get_xunits, set_xunits, 
                      doc = 'x-axis list of units for each dimension')
    yunit = property(get_yunit, set_yunit, 
                     doc = 'y-axis unit')

    def __getitem__(self, index):
        if type(index) in (types.IntType, slice):
            if type(index) is types.IntType and len(self.shape) == 1:
                return self._y[index]
            else:
                index = (index,)

        ## Extend index from left to have same length as dimension
        index = (self.ndim - len(index)) * [slice(None)] + list(index)

        if len(index) > self.ndim:
            raise IndexError('Index order exceeds the number of dimensions')

        newxlist = [x[indexpart] 
                    for x, indexpart in zip(self._xlist, index)
                    if type(indexpart) is not types.IntType]
        newxlabels = [label for label,indexpart in zip(self.xlabels, index)
                      if type(indexpart) is not types.IntType]
        newxunits = [label for label,indexpart in zip(self.xunits, index)
                     if type(indexpart) is not types.IntType]

        newy = self._y[index]


        return Waveform(tuple(newxlist), newy, 
                        xlabels = newxlabels,
                        xunits = newxunits,
                        ylabel = self.ylabel,
                        yunit = self.yunit)

    def getitem(self, key):
        """Get item of each Y-value

        >>> wave = Waveform([[1,2]], array([{'a':1, 'b':2}, {'a':3, 'b':4}]))
        >>> wave.getitem('a')
        Waveform(array([1, 2]), array([1, 3]))

        """
        def getitem(d):
            return d[key]
        
        ufuncgetitem = np.vectorize(getitem)
        return Waveform(self._xlist, ufuncgetitem(self._y), 
                        xlabels = self.xlabels)

    def __repr__(self):
        if self._dim > 1:
            xlist = self._xlist
        else:
            xlist = self._xlist[0]
        return self.__class__.__name__ + "(" + repr(xlist) + ", " + \
            repr(self.y) + ")"

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
                            

## Utility functions
def iswave(w):
    """Returns true if argument is a waveform"""
    return isinstance(w, Waveform)

def assert_waveform(w):
    assert iswave(w), "%s is not a waveform object"%str(w)

def applyfunc(func, w, funcname = None):
    if iswave(w):
        outw  = Waveform(w.x, w.y, xlabels=w.xlabels, xunits=w.xunits)
        outw.y = func(outw._y)
        if w.ylabel:
            if funcname:
                outw.ylabel = funcname + '(' + w.ylabel + ')'
            else:
                outw.ylabel = func.__name__ + '(' + w.ylabel + ')'
        return outw
    else:
        return func(w)

def applyfunc_and_reducedim(func, w, axis = -1, ylabel = None, yunit = None):
    """Apply a function that reduces the dimension by one and return a new waveform or float if zero-rank
    
    """
    newyshape = list(w._y.shape)
    del newyshape[axis]
    newy = apply_along_axis(func, axis, w._y).reshape(newyshape)

    if ylabel != None:
        ylabel = func.__name__ + '(' + ylabel + ')'
    
    return reducedim(w, newy, axis=axis, ylabel=ylabel, yunit=yunit)

def reducedim(w, newy, axis=-1, ylabel=None, yunit=None):
    """Reduce the dimension by one and return a new waveform or float if zero-rank"""

    if rank(newy) == 0:
        return np.asscalar(newy)

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

    Examples:
    
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
    from pycircuit.utilities import rst

    xvalues = cartesian(waveforms[0]._xlist)
    yvalues = zip(*[list(w._y.flat) for w in waveforms])

    xlabels = waveforms[0].xlabels
    ylabels = [w.ylabel for w in waveforms]
    xunits = waveforms[0].xunits
    yunits = [w.yunit for w in waveforms]

    hasunits = not reduce(operator.__and__, [yunit == '' for yunit in yunits])
    
    if hasunits:
        return rst.table(map(lambda x,y: list(x) + list(y), 
                                       [xlabels] + [xunits] + xvalues, 
                                       [ylabels] + [yunits] + yvalues), 
                              headerrows = 2)
    else:
        return rst.table(map(lambda x,y: list(x) + list(y), 
                                       [xlabels] + xvalues, 
                                       [ylabels] + yvalues))
    

def apply_along_axis_with_idx(func1d,axis,arr,*args):
    """ Execute func1d(arr[i], i, *args) where func1d takes 1-D arrays
        and arr is an N-d array.  i varies so as to apply the function
        along the given axis for each 1-d subarray in arr.
    """
    arr = np.asarray(arr)
    nd = arr.ndim
    if axis < 0:
        axis += nd
    if (axis >= nd):
        raise ValueError("axis must be less than arr.ndim; axis=%d, rank=%d."
            % (axis,nd))
    ind = [0]*(nd-1)
    i = np.zeros(nd,'O')
    indlist = range(nd)
    indlist.remove(axis)
    i[axis] = slice(None,None)
    outshape = np.asarray(arr.shape).take(indlist)
    i.put(indlist, ind)
    res = func1d(arr[tuple(i.tolist())], tuple(i.tolist()), *args)
    #  if res is a number, then we have a smaller output array
    if isscalar(res):
        outarr = np.zeros(outshape,np.asarray(res).dtype)
        outarr[tuple(ind)] = res
        Ntot = np.product(outshape)
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
        Ntot = np.product(outshape)
        holdshape = outshape
        outshape = list(arr.shape)
        outshape[axis] = len(res)
        outarr = np.zeros(outshape,np.asarray(res).dtype)
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

    Examples:

    >>> w1 = Waveform(array([1,2,3]),array([3,5,6]))
    >>> w2 = Waveform(array([1,2,3]),array([1,1.5,2]))
    >>> w3 = Waveform(array([1,2]),array([3,5,6]))
    >>> compatible(w1, w2)
    True
    >>> compatible(w1, w3)
    False
    
    """
    return set([w.shape for w in args]) == set([args[0].shape])

def compose(wlist, x = None, xlabel = None):
    """Compose list of waveforms into a new waveform where the 
    waveform list becomes the outer sweep

    Examples:

    >>> wlist=[Waveform(array([1,2,3]),array([3,5,6])), \
               Waveform(array([1,2,3]),array([1,1.5,2]))]
    >>> w = compose(wlist, x = array([1,2]), xlabel = 'index')
    >>> w
    Waveform([array([1, 2]), array([1, 2, 3])], array([[ 3. ,  5. ,  6. ],
           [ 1. ,  1.5,  2. ]]))
    >>> w.xlabels
    ('index', 'x0')
    
    """
    if not compatible(*wlist):
        return ValueError('Waveforms in wlist are not compatible')
    if x != None and len(wlist) != len(x):
        return ValueError('Number of x-values must be the same '
                          'as the number of waveforms')

    newy = np.array([w.y for w in wlist])

    if x == None:
        newx = [range(len(wlist))] + wlist[0].x
    else:
        newx = [x] + wlist[0].x

    if xlabel == None:
        xlabel = 'composite index'

    return Waveform(newx, newy,
                    xlabels = [xlabel] + list(wlist[0].xlabels), 
                    ylabel = wlist[0].ylabel,
                    xunits = [''] + list(wlist[0].xunits), 
                    yunit = wlist[0].yunit)

# Cartesian product operator of list of lists
def cartesian(listList):
    if listList:
        result = []
        prod = cartesian(listList[:-1])
        for x in prod:
            for y in listList[-1]:
                result.append(x + (y,))
        return result
    return [()]

def onedim_index(index, axis, ndim):
    """Return an index from a 1-dimensional index along the given axis
    
    >>> index = onedim_index(np.index_exp[1,2], 1, 3)
    >>> A=np.eye(3)
    >>> A[index]
    array([[ 0.,  0.],
           [ 1.,  0.],
           [ 0.,  1.]])
    
    """
    
    return (slice(None),) * (axis % ndim) + (index,)

            

def wavefunc(func):
    """Decorator for creating free functions from waveform methods
    
    If the first argument is a waveform the waveform method with the same
    name as the function will be called otherwise the decorated function 
    is called instead.

    """
    def g(*args, **kvargs):
        if iswave(args[0]):
            return getattr(Waveform, func.__name__)(*args, **kvargs)
        else:
            return func(*args, **kvargs)

    g.__name__ = func.__name__

    ## Copy docstring from waveform method
    if func.__doc__:
        g.__doc__ = func.__doc__
    else:
        g.__doc__ = getattr(Waveform, func.__name__).__doc__

    return g

if __name__ == "__main__":
    import doctest
    doctest.testmod()
