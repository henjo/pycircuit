import result

class InternalResultSet(result.ResultSet):
    """ResultSet implementation where the results are stored in memory in a dictionary

    >>> rs = InternalResultSet()
    >>> rs['test'] = InternalResult(signals = {'testA': 1})
    >>> print rs
    
    """
    def __init__(self, results = {}):
        self.results = dict(results)

    def __setitem__(self, name, result):
        self.results[name] = result
        
    def ___len__(self):
        return len(self.results)

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
    
    def __len__(self):
        return len(self.signals)

    def __setitem__(self, name, signal):
        self.results[name] = signal

    def getSignalNames(self): return self.signals.keys()
    def getSignal(self, name): return self.signals[name]
    
    def storeSignal(self, name, signal):
        self.signals[name] = signal
        self.o.update()

if __name__ == "__main__":
    import doctest
    doctest.testmod()
