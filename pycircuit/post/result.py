# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import pycircuit.utilities.rst as rst

class ResultDict(object):
    """The ResultDict class manages results from a circuit simulation.

    """

    def __dict__(self):
        return dict([(k, self[k]) for k in self.keys()])
    def __delitem__(self, key):
        raise NotImplementedError()
    def __getitem__(self, key):
        raise NotImplementedError()
    def __setitem__(self, key, value):
        raise NotImplementedError()
    def __len__(self):
        raise NotImplementedError()
    def __contains__(self, key):
        return key in self.keys()
    def keys(self):
        """List of ResultDict's keys"""
        raise NotImplementedError()

    def values(self):
        return [self[k] for k in self.keys()]

    def __str__(self):
        s = 'ResultDict: %s'%(str(self.keys()))
        return s
    # + rsttable.toRSTtable([['Result name', '# of signals']] + [[resultname, str(len(self[resultname]))] for resultname in self.keys()])
    def __repr__(self):
        return self.__class__.__name__ + '(' + str(list(self.keys())) + ')'

    @property
    def table(self):
        return rsttable.toRSTtable([['Signal', 'Value']] + [(k, str(v)) for k,v in dict(self).items()])

class IVResultDict(ResultDict):
    """Result dictionary for storing voltages and currents"""
    
    def v(self, plus, minus=None):
        """Returns the voltage between the plus and minus node or potential of plus node"""
        raise NotImplementedError()

    def i(self, terminal):
        """Returns the current flowing into the given terminal or branch object"""
        raise NotImplementedError()
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
