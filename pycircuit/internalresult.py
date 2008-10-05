import result

class InternalResultSet(result.ResultSet):
    """ResultSet implementation where the results are stored in memory in a dictionary

    >>> rs = InternalResultSet()
    >>> rs.storeResult('test', InternalResult(signals = {'testA': 1}))
    >>> rs
    
    """
    def __init__(self, results = {}):
        self.results = dict(results)
        
    def getResultNames(self):
        return self.results.keys()

    def getResult(self, name):
        return self.results[name]

    def storeResult(self, name, result):
        self.results[name] = result

class InternalResult(result.Result):
    def __init__(self, signals = {}):
        self.signals = dict(signals)
        result.Result.__init__(self)
    
    def getSignalNames(self): return self.signals.keys()
    def getSignal(self, name): return self.signals[name]
    
    def storeSignal(self, name, signal):
        self.signals[name] = signal
        self.o.update()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
