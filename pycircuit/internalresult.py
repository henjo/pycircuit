import result

class InternalResultSet(result.ResultSet):
    """ResultSet implementation where the results are stored in memory in a dictionary"""
    def __init__(self, results = {}):
        self.results = results
        
    def getResultNames(self):
        return self.results.keys()

    def getResult(self, name):
        return self.results[name]

    def storeResult(self, name, result):
        self.results[name] = result

class InternalResult(result.Result):
    def __init__(self, signals = {}):
        self.signals = signals
        result.Result.__init__(self)
    
    def getSignalNames(self): return self.signals.keys()
    def getSignal(self, name): return self.signals[name]
    
    def storeSignal(self, name, signal):
        self.signals[name] = signal
        self.o.update()

