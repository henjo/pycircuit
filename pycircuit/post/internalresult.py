# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import result

class InternalResultDict(result.ResultDict):
    """ResultSet implementation where the results are stored in memory in a dictionary

    >>> rs = InternalResultDict()
    >>> rs['test'] = InternalResultDict({'testA': 1})
    >>> print rs
    ResultDict: ['test']
    
    """
    def __init__(self, items = {}):
        self.items = dict(items)

    def keys(self): return self.items.keys()

    def __delitem__(self, key):
        del self.items[key]

    def __getitem__(self, key):
        return self.items[key]

    def __setitem__(self, key, value):
        self.items[key] = value
        
    def __len__(self):
        return len(self.results)


if __name__ == "__main__":
    import doctest
    doctest.testmod()
