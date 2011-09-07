# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from numpy.testing import assert_equal, assert_array_almost_equal, assert_array_equal
import unittest

from pycircuit.utilities.misc import *

def test_create_index_vectors():
    indices = [0, 1, 1, 2, 2, 2]
    src = N.array([1,2,3,4,5,6])

    refdest = N.zeros(3)
    inplace_add_selected_ref(refdest, indices, src)

    dest = N.zeros(3)

    for dst_i, src_i in create_index_vectors(indices): 
        dest[dst_i] += src[src_i]

    assert_array_equal(refdest, dest)

def test_inplace_add_selected():
    src = N.array([1,2,3,9,-1,1])
    indices = [0, 1, 1, 2, 2, 2]

    indices_list = create_index_vectors(indices)
    
    a = N.zeros(3)
    inplace_add_selected(a, indices_list, src)

    aref = N.zeros(3)
    inplace_add_selected_ref(aref, indices, src)

    assert_array_equal(a, aref)

@unittest.skip("Skip failing test")
def test_ObserverSubject_init():
    a = ObserverSubject()
    assert_equal(a._observers,[])

@unittest.skip("Skip failing test")
def test_ObserverSubject_notify():
    
    class MySubject(ObserverSubject):
        def hitme(self):
            self.notify(args=10)

    class objectWithUpdateMethod(object):
        def __init__(self):
            self.reset()

        def reset(self):
            self.value = None
            self.value1 = None
            
        def update(self, subject, args):
            self.value = ('update', subject, args)

        def update1(self, subject, args):
            self.value1 = ('update1', subject, args)
            
    a = MySubject()
    b = objectWithUpdateMethod()
    
    a.attach(b)
    a.attach(b, 'update1')

    a.hitme()
    
    assert_equal(b.value, ('update', a, 10))
    assert_equal(b.value1, ('update1', a, 10))

    b.reset()

    # Test detach
    a.detach(b, 'update')
    a.hitme()

    assert_equal(b.value, None)
    assert_equal(b.value1, ('update1', a, 10))

    a.detach(b, 'update1')
    a.hitme()

    assert_equal(b.value, None)
    assert_equal(b.value1, ('update1', a, 10))
