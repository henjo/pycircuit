import rsttable

class ResultSet(object):
    """The ResultCollection class handles a set of results"""

    def __getitem__(self, name):
        return self.getResult(name)
    def __setitem__(self, name, result):
        raise Exception('Not implemented')
    def ___len__(self):
        raise Exception('Not implemented')
    def keys(self):
        return self.getResultNames()

    def getResultNames(self):
        pass
    def getResult(self, name):
        pass

    def __str__(self):
        s = 'Resultset results:\n'
        return rsttable.toRSTtable([['Result name', '# of signals']] + [[resultname, str(len(self[resultname]))] for resultname in self.keys()])
    
        return self.__class__.__name__ + '(' + str(list(self.getResultNames())) + ')'


class Result(object):
    """The result class manages results from a circuit simulation.

       A result contains the values of one or more signals. Result can
       handle multi dimensional sweeps but not several sweeps.
    """
    def __getitem__(self, name):
        return self.getSignal(name)
    def __setitem__(self, name, result):
        raise Exception('Not implemented')
    def keys(self):
        return self.getSignalNames()
    
    def getSweepValues(self, dimension):
        """Get a numpy array of sweep values from sweep dimension dimension"""
        raise Exception('Not implemented')

    def getSignalNames(self):
        """Returns a tuple of available signals in the result"""
        raise Exception('Not implemented')

    def ___len__(self):
        raise Exception('Not implemented')
    
    def getSignal(self, name):
        """Returns a Waveform object if the parameter is swept and a scalar otherwise.
           If name is None it will return a dictionary of all signals keyed by
           the signal names
        """
        raise Exception('Not implemented')

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
