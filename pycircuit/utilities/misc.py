import numpy as N
from itertools import izip, groupby

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
    """
    Adds values to selected elements of dest

    The difference compared to dest[indices] += values is that indices may contain
    duplicate indices which would add that element several times.

    >>> a = N.zeros(3)
    >>> inplace_add_selected(a, [0,1,1], N.array([1,2,3]))
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
    """Adds values to selected elements of dest using precalculated indices lists
    
    The indiceslist is created as:
    
    indiceslist = create_index_vectors(indices)
    
    """
    for dst_i, src_i in indiceslist: 
        for dst_j, src_j in indiceslist: 
            dest[[[i] for i in dst_i], dst_j] += values[[[i] for i in src_i], src_j]

def create_index_vectors(indices):
    """Create list of index vectors suitable for use repeated calls to in place addition
    
    >>> indices = [0, 1, 1, 2, 2, 2]
    >>> dest = N.zeros(3)
    >>> src = N.array([1,2,3,4,5,6])
    >>> inplace_add_selected_ref(dest, indices, src)
    >>> ref = dest
    >>> dest = N.zeros(3)
    >>> for dst_i, src_i in create_index_vectors(indices): dest[dst_i] += src[src_i]
    >>> N.alltrue(ref == dest)
    True
    
    """

    i = 0
    dupes = {}
    for k,group in groupby(indices):
        n = len(list(group))

        dupes[k] = i, n

        i += n
        
    result = []

    while len(dupes) > 0:
        src_i = []
        dst_i = []
        for k, v in dupes.items():
            i, n = v
            
            if n > 1:
                dupes[k] = i+1, n-1
            else:
                del dupes[k]

            src_i.append(i)
            dst_i.append(k)

        result.append((dst_i, src_i))
        


    return result
    
    
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
