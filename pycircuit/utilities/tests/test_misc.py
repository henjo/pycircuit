# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from numpy.testing import assert_equal, assert_array_almost_equal, assert_array_equal

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
        
def test_ObserverSubject_init():
    a = ObserverSubject()
    assert_equal(a._observers,[])

def test_ObserverSubject_attach():
    a = ObserverSubject()
    a.attach('APA')
    assert_equal(a._observers,['APA'])

def test_ObserverSubject_detach():
    a = ObserverSubject()
    a._observers = ['APA']
    assert_equal(a._observers,['APA'])    
    a.detach('APA')
    assert_equal(a._observers,[])

def test_ObserverSubject_attach_detach():
    a = ObserverSubject()
    a.attach('APA')
    assert_equal(a._observers,['APA'])
    a.detach('APA')
    assert_equal(a._observers,[])

def test_ObserverSubject_notify():
    
    class objectWithUpdateMethod(object):
        def __init__(self):
            self.value = 'Not updated'

        def update(self,input):
            self.value = 'Updated'
            
    a = ObserverSubject()
    b = objectWithUpdateMethod()
    
    a.attach(b)
    a.notify()
    
    assert_equal(b.value,'Updated')
