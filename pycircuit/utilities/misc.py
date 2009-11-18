# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import numpy as N
from itertools import izip, groupby
from operator import itemgetter
import tempfile
import shutil

def indent(s, n=4, notfirstline = False):
    """Indent string

    >>> indent("apa\\nrapa\\nbapa", 4)
    '    apa\\n    rapa\\n    bapa'
    >>> indent("apa\\nrapa\\nbapa", 4, notfirstline=True)
    'apa\\n    rapa\\n    bapa'

    """
    if notfirstline:
        return ('\n' + n*' ').join(s.split('\n'))
    else:
        return '\n'.join([n*' ' + line for line in s.split('\n')])

def inplace_add_selected_ref(dest, indices, values):
    """Adds an array to the selected indices of a destination array

    The difference compared to dest[indices] += values is that indices may contain
    duplicate indices which would add that element several times.

    >>> a = N.zeros(3)
    >>> inplace_add_selected_ref(a, [0,1,1], N.array([1,2,3]))
    >>> a
    array([ 1.,  5.,  0.])

    """

    for i, v in izip(indices, values):
        dest[i] += v
        

def inplace_add_selected(dest, indiceslist, values):
    """Adds values to selected elements of dest using precalculated indices lists
    
    The indiceslist is created as:
    
    indiceslist = create_index_vectors(indices)
    
    """
    for dst_i, src_i in indiceslist: 
        dest[dst_i] += values[src_i]

def inplace_add_selected_2d(dest, indiceslist, values):
    """Adds a 2-dimensional array to the selected indices of a 2-d destination array
    
    The indiceslist is created as:
    
    indiceslist = create_index_vectors(indices)
    
    """
    for dst_i, src_i in indiceslist: 
        for dst_j, src_j in indiceslist: 
            dest[[[i] for i in dst_i], dst_j] += values[[[i] for i in src_i], src_j]

def create_index_vectors(indices):
    """Create list of index vectors suitable for use repeated calls to in place addition
    
    >>> create_index_vectors([2, 0, 1, 0])
    [([0, 1, 2], [3, 2, 0]), ([0], [1])]
    
    """

    sorted_indices = sorted(enumerate(indices), key=itemgetter(1))

    dupes = {}
    for dst_i, group in groupby(sorted_indices, key=itemgetter(1)):
        src_indices = [src_i[0] for src_i in group]
        dupes[dst_i] = src_indices

    result = []
    while len(dupes) > 0:
        src_i_list = []
        dst_i_list = []
        for dst_i, src_indices in dupes.items():
            src_i = src_indices.pop()

            if len(src_indices) == 0:
                del dupes[dst_i]

            src_i_list.append(src_i)
            dst_i_list.append(dst_i)

        result.append((dst_i_list, src_i_list))

    return result
    
def combinations(iterable, r):
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = list(range(r))
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)

class ObserverSubject(object):
    """Subject class in the observer design pattern"""
    def __init__(self):
        self._observers = []

    def attach(self, observer):
        """Attach observer that will be notified when an item changes"""
        if not observer in self._observers:
            self._observers.append(observer)

    def detach(self, observer):
        try:
            self._observer.remove(observer)
        except ValueError:
            pass

    def notify(self, modifier=None):
        for observer in self._observers:
            if modifier != observer:
                observer.update(self)

class TempDir(object):
    def __init__(self, srcdir=None, keep=False):
        self.keep = keep
        self.dir = tempfile.mkdtemp()

        if srcdir:
            copytree(srcdir, self.dir)

    def __str__(self):
        return self.dir
    def __del__(self):
        if not self.keep:
            shutil.rmtree(self.dir)



if __name__ == "__main__":
    import doctest
    doctest.testmod()
