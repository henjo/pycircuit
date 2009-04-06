# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from numpy.testing import assert_array_almost_equal, assert_array_equal

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
        
    


    
