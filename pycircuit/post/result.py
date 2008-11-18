import pycircuit.utilities.rst as rst

class ResultDict(object):
    """The ResultDict class manages results from a circuit simulation.

    """

    def __dict__(self):
        return dict([(k, self[k]) for k in self.keys()])
    def __getitem__(self, key):
        raise Exception('Not implemented')
    def __setitem__(self, key, value):
        raise Exception('Not implemented')
    def __len__(self):
        raise Exception('Not implemented')
    def __contains__(self, key):
        return key in self.keys()
    def keys(self):
        raise Exception('Not implemented')

    def __str__(self):
        s = 'ResultDict: %s'%(str(self.keys()))
        return s
    # + rsttable.toRSTtable([['Result name', '# of signals']] + [[resultname, str(len(self[resultname]))] for resultname in self.keys()])
    def __repr__(self):
        return self.__class__.__name__ + '(' + str(list(self.keys())) + ')'


    @property
    def table(self):
        return rsttable.toRSTtable([['Signal', 'Value']] + [(k, str(v)) for k,v in dict(self).items()])
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
