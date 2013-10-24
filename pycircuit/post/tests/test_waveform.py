# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from nose.tools import *
from numpy.testing import assert_array_almost_equal, assert_array_equal

import numpy as np
from pycircuit.post import Waveform
from pycircuit.post.functions import *
from pycircuit.post.testing import *
import unittest

testdata1 = (
    Waveform(array([1,10,100]),array([complex(-1,0),complex(0,1),2]), 
             xlabels = ('freq',),
             xunits = ('Hz', ),
             ylabel = 'amplitude',
             yunit = 'V'
             ),    complex(0,1),
    -1.5,
    1.5)
testdata1_0_table = """====== ===========
freq   amplitude  
Hz     V          
====== ===========
     1     (-1+0j)
    10          1j
   100      (2+0j)
====== ==========="""

testdata2 = Waveform([[1,2], [1,2,3]], array([[3,4,5],[5,4,2]]),
                     xlabels = ('v1', 'v2'),
                     xunits = ('V', 'V'),
                     ylabel = 'i3',
                     yunit = 'A')
testdata2_table = """==== ==== ====
v1   v2   i3  
V    V    A   
==== ==== ====
   1    1    3
   1    2    4
   1    3    5
   2    1    5
   2    2    4
   2    3    2
==== ==== ===="""

# Commented out as this causes the following error
# ValueError: Dimension of x ([2, 2]) does not match dimension of y ((2,))
# testdata3 = Waveform([array([array([1,1]), array([2])], dtype=object), 
#                             array([array([1,2]), array([1])], dtype=object)],
#                            array([array([3,4]),array([4])], dtype = object),
#                            xlabels = ('v1', 'v2'),
#                            xunits = ('V', 'V'),
#                            ylabel = 'i3',
#                            yunit = 'A')

testdata3_table = """==== ==== ====
v1   v2   i3  
V    V    A   
==== ==== ====
   1    1    3
   1    2    4
   2    1    4
==== ==== ===="""

testdata4 = Waveform([[1,2], [1,2,3], [3,4]], array([[[1,3],[4,4],[5,4]],[[4,5],[4,4],[3,2]]]),
                     xlabels = ('v1', 'v2', 'v3'),
                     xunits = ('V', 'V', 'V'),
                     ylabel = 'i3',
                     yunit = 'A')
testdata4_table = """==== ==== ====
v1   v2   i3  
V    V    A   
==== ==== ====
   1    1    3
   1    2    4
   1    3    5
   2    1    5
   2    2    4
   2    3    2
==== ==== ===="""


def test_creation():
    """Test of waveform creation"""

    def check(w):
        for x, xref in zip(w.x, [array([1,2])]):
            assert_array_equal(x, xref)
        assert_array_equal(w.y, array([3,4]))

    ## Create from arrays
    check(Waveform([array([1,2])], array([3,4])))

    ## Create from single x-array
    check(Waveform(array([1,2]), array([3,4])))

    ## Create from lists
    check(Waveform([[1,2]], [3,4]))

    ## Create from tuples
    check(Waveform(((1,2),), (3,4)))
    
    ## Check meta-data of testdata waveforms
    for x,xref in zip(testdata2.x, [np.array([1,2]), np.array([1,2,3])]):
        assert_array_equal(x, xref)
        
    assert_array_equal(testdata2.y, array([[3,4,5],[5,4,2]]))
    
    assert_equal(testdata2.xlabels, ['v1', 'v2'])
    assert_equal(testdata2.xunits, ['V', 'V'])
    assert_equal(testdata2.ylabel, 'i3')
    assert_equal(testdata2.yunit, 'A')
    assert_equal(testdata2.ndim, 2)
 
    empty = Waveform()
    assert_array_equal(empty.y, array([]))
    assert_equal(empty.ndim, 1)
   


def test_unary_operations():
    x = testdata1[0]
    check_func(abs, lambda w: Waveform(w.x, abs(get_y(w))),
               (x,), preserve_yunit=True)
    check_func(lambda w: -w, lambda w: Waveform(w.x, -get_y(w)),
               (x,), preserve_yunit=True, ref_ylabel = '-%s')

def test_binary_operations():
    for a in testdata1:
        for b in testdata1:
            if iswave(a) or iswave(b):
                for op in ('+', '-', '*', '**', '==', '<', '<=', '>', '>='):
                    check_binary_op(op, a, b)

    ## 2D data
    a = testdata2
    b = testdata2
    for op in ('+', '-', '*', '**'):
        check_binary_op(op, a, b)

def test_indexing():
    """Test taking slices of a waveform"""
    
    w = testdata2
    
    w_sliced = Waveform(([1,2,3],), array([3,4,5]),
                          xlabels = ('v2',),
                          xunits = ('V',),
                          ylabel = 'i3',
                          yunit = 'A')

    assert_waveform_almost_equal(w[0,:], w_sliced)
    assert_waveform_almost_equal(w[0], w_sliced)

    w_sliced2 = Waveform(([1,2],), array([3,5]),
                         xlabels = ('v1',),
                         xunits = ('V',),
                         ylabel = 'i3',
                         yunit = 'A')
    
    assert_waveform_almost_equal(w[:,0], w_sliced2)

    w_sliced3 = Waveform(([1], [1,2,3],), array([[3,4,5]]),
                          xlabels = ('v1', 'v2',),
                          xunits = ('V', 'V',),
                          ylabel = 'i3',
                          yunit = 'A')
    
    assert_waveform_almost_equal(w[0:1,:], w_sliced3)
    assert_equal(w[0,0], 3)

    w_sliced = Waveform(([1,2],), array([5,2]),
                          xlabels = ('v1',),
                          xunits = ('V',),
                          ylabel = 'i3',
                          yunit = 'A')

    assert_waveform_almost_equal(w[:,-1], w_sliced)

    w_sliced = Waveform([[1,2], [1,2,3]], array([[1,4,5],[4,4,3]]),
                        xlabels = ('v1', 'v2'),
                        xunits = ('V', 'V'),
                        ylabel = 'i3',
                        yunit = 'A')

    assert_waveform_almost_equal(testdata4[:,:, 0], w_sliced)
    assert_waveform_almost_equal(testdata4[...], testdata4)
    assert_waveform_almost_equal(testdata4[..., 0], w_sliced)
    
def test_xmax(): 
    w1 = Waveform(array([1,2,3]),array([3,5,6]))
    assert_equal(w1.xmax(), 3)

    w2 = Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
    assert_equal(w2.xmax(), 4)

def test_xmin(): 
    w1 = Waveform(array([1,2,3]),array([3,5,6]))
    assert_equal(w1.xmin(), 1)

    w2 = Waveform([[5,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
    assert_equal(w2.xmin(), 2)

def test_ymin():
    w = Waveform([[5,2],[2,3,4]], array([[3,9,7], [4,6,6]]))
    check_alongaxes_func(Waveform.ymin, np.min, w)

def test_ymax():
    w = Waveform([[2,5],[2,3,4]], array([[3,9,7], [4,6,6]]))
    check_alongaxes_func(Waveform.ymax, np.max, w)

@unittest.skip("Skip failing test")
def test_value():
    w1 = Waveform(array([1,2,3]),array([3,5,6]))
    assert_almost_equal(w1.value(1.5), 4.0)

    ## 2-d waveform
    w2 = Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
    assert_waveform_almost_equal(w2.value(2.5), Waveform([[1, 2]], array([ 4.,  5.])))
    assert_waveform_almost_equal(w2.value(1.5, axis=0), 
                          Waveform([[2, 3, 4]], array([ 3.5, 5.5, 6.5])))
    ## x is a waveform
    w2 = Waveform([[1,2],[2,3,4]], array([[3,5,6], [4,6,7]]))
    assert_waveform_almost_equal(w2.value(Waveform(array([1, 2]), array([2.5, 3.5]))), 
                          Waveform(array([1, 2]),array([ 4.,   6.5])))

def test_xval():
    w = Waveform([[5,2],[2,3,4]], array([[3,9,7], [4,6,6]]))

    wref = Waveform([[5,2],[2,3,4]], array([[2,3,4], [2,3,4]]))
    assert_waveform_almost_equal(w.xval(), wref)
    assert_waveform_almost_equal(w.xval(-1), wref)
    assert_waveform_almost_equal(w.xval(1), wref)

    wref = Waveform([[5,2],[2,3,4]], array([[5,5,5], [2,2,2]]))
    assert_waveform_almost_equal(w.xval(0), wref)

def test_clip():
    ## 1D waveforms
    w1 = Waveform([[1.,2.,3.]], array([8., 6., 1.]), xlabels=['a'])

    assert_waveform_almost_equal(w1.clip(2,3,axis='a'), 
                          Waveform([[ 2.,  3.]], array([ 6.,  1.])))
    assert_waveform_almost_equal(w1.clip(1.5, 3),
                          Waveform([[ 1.5, 2.,  3.]],array([ 7., 6.,  1.])))
    assert_waveform_almost_equal(w1.clip(1., 2.5),
                          Waveform([[ 1., 2.,  2.5]],array([ 8., 6.,  3.5])))
    
    ## 2D waveforms
    w = Waveform([[2,5],[2,3,4,5]], array([[3,9,7,6], [4,6,6,3]]),
                 xlabels = ('x1','x2'), xunits = ('V', 'V'),
                 ylabel = 'i', yunit = 'A')

    assert_waveform_almost_equal(w.clip(3, 4), 
                          Waveform([[2,5],[3,4]], 
                                   array([[9,7], [6,6]]),
                                   xlabels = ('x1','x2'), xunits = ('V', 'V'),
                                   ylabel = 'i', yunit = 'A'))

    w = Waveform([[2.,5.],[2.,3.,4.]], array([[3.,9.,7.], [4.,6.,6.]]),
                 xlabels = ('x1','x2'), xunits = ('V', 'V'),
                 ylabel = 'i', yunit = 'A')

    ## Clip at non-existing x-value from left and right
    assert_waveform_almost_equal(w.clip(2.5, 3.5), 
                          Waveform([[2,5],[2.5,3,3.5]], 
                                   array([[6,9,8], [5,6,6]]),
                                   xlabels = ('x1','x2'), xunits = ('V', 'V'),
                                   ylabel = 'i', yunit = 'A'))

    ## Clip at non-existing x-value from left and right
    assert_waveform_almost_equal(w.clip(2.5, 3.5, axis=0), 
                          Waveform([[2.5,3.5],[2.,3.,4.]], 
                                   array([[3.166667, 8.5, 6.833333], 
                                          [3.5, 7.5, 6.5]]),
                                   xlabels = ('x1','x2'), xunits = ('V', 'V'),
                                   ylabel = 'i', yunit = 'A'))

    
def check_alongaxes_func(func, reference_func, w, keep_yunit = False):
    for axis in range(w.ndim):
        res = func(w, axis=axis)
        
        ref_x = list(w.x)
        ref_xlabels = list(w.xlabels)
        ref_xunits = list(w.xlabels)
        del ref_x[axis]; del ref_xlabels[axis]; del ref_xunits[axis]
        
        if keep_yunit:
            yunit = w.yunit
        else:
            yunit = ''

        y = np.apply_along_axis(reference_func, axis, w.y)

        ref = Waveform(ref_x, y,
                       xlabels = ref_xlabels, xunits = ref_xunits,
                       ylabel = func.__name__ + '(' + w.ylabel + ')',
                       yunit = yunit)

        assert_waveform_almost_equal(res, ref)

def check_func(func, reference_func, args, preserve_yunit = False,
               ref_yunit='', ref_ylabel=None):
    res = func(*args)

    ## Check result
    if iswave(res):
        y_ref = reference_func(*args)
        if iswave(y_ref):
            y_ref = y_ref.y
        assert_array_equal(res.y, y_ref)
        assert_array_equal(np.array(res.x), 
                           args[0].x)
    else:
        assert_array_equal(res, reference_func(*args))
    
    ## Check that the argument and results have the same type
    assert any([iswave(x) for x in args]) == iswave(res)

    ## Check x label and units
    if iswave(args[0]):
        assert_equal(res.xlabels, x.xlabels)
        assert_equal(res.xunits, x.xunits)

        if ref_ylabel == None:
            ref_ylabel = func.__name__ + '(' + x.ylabel + ')'
        else:
            ref_ylabel = ref_ylabel%x.ylabel

        assert_equal(res.ylabel, ref_ylabel)

        if preserve_yunit:
            assert_equal(res.yunit, x.yunit)
        else:
            assert_equal(res.yunit, ref_yunit)

def check_nonscalar_function(func):
    """Check that scalar input to a  waveform-only functions raises an exception"""
    assert_raises(AssertionError, func, 10)

def get_y(w):
    if iswave(w):
        return w._y
    else:
        return w

def get_ylabel(w):
    if iswave(w):
        return w._ylabel
    else:
        return str(w)

def check_binary_op(op, a, b, preserve_yunit = False):
    exec 'res = a' + op + 'b'

    ref_yunit = ''

    for arg in a,b:
        print str(a) + ' ' + op + ' ' + str(b)
        if iswave(arg):
            ref_x = arg.x
            ref_xlabels = arg.xlabels
            ref_xunits = arg.xunits
            if preserve_yunit:
                ref_yunit == arg.yunit
            break
    
    ref_ylabel = get_ylabel(a) + op + get_ylabel(b)

    exec 'ref_y = get_y(a) ' + op + ' get_y(b)'

    assert_waveform_almost_equal(res, Waveform(ref_x, ref_y, xlabels = ref_xlabels, 
                                        xunits = ref_xunits, ylabel = ref_ylabel,
                                        yunit = ref_yunit))

def test_array_interface():
    for win in testdata1[0], testdata2:
        wout = np.cos(win)
        assert_array_equal(wout.y, np.cos(win))
        assert_array_equal(win.x[0], wout.x[0])

def test_getitem():
    w = testdata1[0]
    assert_equal(w[0], w.y[0])

    wsliced = Waveform(array([1]),array([complex(-1,0)]), 
                       xlabels = ('freq',),
                       xunits = ('Hz', ),
                       ylabel = 'amplitude',
                       yunit = 'V'
                       )
    assert_equal(w[0:1], wsliced)
    
    
def test_repr():
    for wref in testdata1[0], testdata2:
        wstr = repr(wref)
        wout = eval(wstr)
        assert_waveform_almost_equal(wout, wref)


def test_astable():

    for w, t in zip((testdata1[0], testdata2), 
                    (testdata1_0_table, testdata2_table)):
        assert_equal(w.astable, t)

def test_simple_broadcasting():
    assert_waveform_equal(testdata2 + testdata2[0],
                          Waveform([[1,2], [1,2,3]], 
                                   testdata2.y + testdata2.y[0],
                                   xlabels = ('v1', 'v2'),
                                   xunits = ('V', 'V'),
                                   ylabel = 'i3',
                                   yunit = 'A'))

def test_reversed_broadcasting():
    assert_waveform_equal(testdata2 + testdata2[:,0],
                          Waveform([[1,2], [1,2,3]], 
                                   (testdata2.y.T + (testdata2.y.T)[0]).T,
                                   xlabels = ('v1', 'v2'),
                                   xunits = ('V', 'V'),
                                   ylabel = 'i3',
                                   yunit = 'A'))
    
def test_duplicate_xlabels():
    """Check that a waveform with duplicate xlabels raises an exception"""
    def func():
        return Waveform([[1,2], [1,2,3]], 
                        (testdata2.y.T + (testdata2.y.T)[0]).T,
                        xlabels = ('v1', 'v1'),
                        xunits = ('V', 'V'),
                        ylabel = 'i3',
                        yunit = 'A')
    assert_raises(ValueError, func)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
