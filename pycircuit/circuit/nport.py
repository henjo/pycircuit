import numpy as N
from circuit import SubCircuit, gnd, R, VS
from analysis import Analysis, AC
from pycircuit.internalresult import InternalResultSet, InternalResult
from copy import copy

class NPort(object):
    """Class to perform calculations and conversion of n-port circuits

    Attributes
    ----------
    n -- number of ports
    K -- (numpy array) transmission parameter matrix, the size is 2*n by 2*n


    """

    def __init__(self, K):
        assert size(K,0) == size(K,1) and k % 2 == 0
        
        self.K = size(K,0)

        self.n = size(K,0) / 2

    def __mul__(self, a):
        """Cascade of two n-ports"""
        return NPort(dot(self.K, a.k))
        
    def __repr__(self):
        return repr(self.K)

class TwoPortAnalysis(Analysis):
    """Analysis to find the 2-ports parameters of a circuit

    The transmission parameters are found as:

    A = v(inp, inn)/v(outp, outn) | io = 0
    B = v(inp, inn)/i(outp, outn) | vo = 0
    C = i(inp, inn)/v(outp, outn) | io = 0
    D = i(inp, inn)/i(outp, outn) | vo = 0

    >>> c = SubCircuit()
    >>> n1 = c.addNode('net1')
    >>> n2 = c.addNode('net2')
    >>> c['R1'] = R(n1, n2, r=9e3)
    >>> c['R2'] = R(n2, gnd, r=1e3)
    >>> res = TwoPort(c, n1, gnd, n2, gnd).run(freqs = N.array([0]))
    >>> res['mu'].y[0]
    (0.1+0j)
    >>> res['gamma'].y[0]
    (0.000111111111111+0j)
    >>> res['zeta'].y[0]
    (1000+0j)
    >>> res['beta'].y[0]
    (1+0j)
    
    """
    
    ACAnalysis = AC
    
    def __init__(self, circuit, inp, inn, outp, outn):
        self.c = circuit

        self.ports = inp, inn, outp, outn
        
    def run(self, freqs, **kvargs):
        result = InternalResult()

        abcd = self.solve(freqs, **kvargs)
        result.storeSignal('ABCD', abcd)

        result.storeSignal('mu', 1/abcd[0,0])
        result.storeSignal('gamma', 1/abcd[0,1])
        result.storeSignal('zeta', 1/abcd[1,0])
        result.storeSignal('beta', 1/abcd[1,1])

        self.result = result

        return result

    def solve(self, freqs, refnode = gnd, complexfreq = False):
        inp, inn, outp, outn = self.ports
                
        ## Add voltage source at input port and create
        ## copies with output open and shorted respectively
        circuit_vs_open = copy(self.c)

        circuit_vs_open['VS_TwoPort'] = VS(inp, inn, v=1.0)

        circuit_vs_shorted = copy(circuit_vs_open)

        circuit_vs_shorted['VL_TwoPort'] = VS(outp, outn, v=0.0)

        ## Run AC-analysis on the two circuits
        ac_open = self.ACAnalysis(circuit_vs_open)
        ac_shorted = self.ACAnalysis(circuit_vs_shorted)

        ac_open.run(freqs, refnode = refnode, complexfreq=complexfreq)

        ac_shorted.run(freqs, refnode = refnode, complexfreq=complexfreq)
        
        A = ac_open.v(inp, inn) / ac_open.v(outp, outn)
        B = ac_shorted.v(inp, inn) / ac_shorted.i('VL_TwoPort.plus')
        C = ac_open.i('VS_TwoPort.minus') / ac_open.v(outp, outn)
        D = ac_shorted.i('VS_TwoPort.minus') / ac_shorted.i('VL_TwoPort.plus')

        return N.array([[A,B],[C,D]])

if __name__ == "__main__":
    import doctest
    doctest.testmod()
