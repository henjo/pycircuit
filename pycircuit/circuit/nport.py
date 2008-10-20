import numpy as N
from circuit import SubCircuit, gnd, R, VS, IS
from analysis import Analysis, AC, Noise
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
        assert N.size(K,0) == N.size(K,1) and N.size(K) % 2 == 0
        
        self.K = K

        self.n = N.size(K,0)

    @classmethod
    def fromY(self, Y):
        """Create a NPort from Y-parameters"""
        d = Y[0,0] * Y[1,1] - Y[0,1] * Y[1,0]
        return NPort(N.array([[-Y[1,1] / Y[1,0], -1.0 / Y[1,0]],
                              [-d / Y[1,0], -Y[0,0] / Y[1,0]]]))

    @classmethod
    def fromZ(self, Z):
        """Create a NPort from Z-parameters"""
        d = Z[0,0] * Z[1,1] - Z[0,1] * Z[1,0]
        return NPort(N.array([[-Z[0,0] / Z[1,0], -d / Z[1,0]],
                              [1.0 / Z[1,0], Z[1,1] / Z[1,0]]]))

    def __mul__(self, a):
        """Cascade of two n-ports"""
        return NPort(dot(self.K, a.k))

    def __floordiv__(self, a):
        """Parallel of two n-ports

        >>> import sympy as S
        >>> a,b,c,d = S.symbols('abcd')
        >>> A = TwoPort(N.array([[a,b], [c,d]]))
        >>> print A // A
        array([[a, 0.5*b],
           [-2c, d]], dtype=object)

        """
        return NPort.fromY(self.Y + a.Y)

    def series(self, a):
        """Series connection with another n-port

        >>> import sympy as S
        >>> a,b,c,d = S.symbols('abcd')
        >>> A = TwoPort(N.array([[a,b], [c,d]]))
        >>> A.series(A).K
        [[ a 2*b]
         [ c/2 d]]
        
        """
        return NPort.fromZ(self.Z + a.Z)

    @property
    def Z(self):
        """Return Z-parameter matrix"""
        if self.n != 2:
            raise Exception('Cannot handle %d-port parameters yet'%self.n)
        
        A = self.K
        d = A[0,0] * A[1,1] - A[0,1] * A[1,0]

        return N.array([[A[0,0] / A[1,0], d / A[1,0]],
                        [1.0 / A[1,0], A[1,1] / A[1,0]]])
    

    @property
    def Y(self):
        """Return Y-parameter matrix"""
        if self.n != 2:
            raise Exception('Cannot handle %d-port parameters yet'%self.n)

        A = self.K
        d = A[0,0] * A[1,1] - A[0,1]*A[1,0]

        return N.array([[A[1,1] / A[0,1], -d / A[0,1]],
                        [-1.0 / A[0,1], A[0,0] / A[0,1]]])
    
        
    def __str__(self):
        return repr(self.K)

class TwoPort(NPort):
    n = 2

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
    >>> res = TwoPortAnalysis(c, n1, gnd, n2, gnd).run(freqs = N.array([0]))
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
    NoiseAnalysis = Noise
    
    def __init__(self, circuit, inp, inn, outp, outn, noise = False, noise_outquantity = 'v'):
        self.c = circuit

        self.ports = inp, inn, outp, outn

        self.noise = noise
        self.noise_outquantity = noise_outquantity
        
    def run(self, freqs, **kvargs):
        result = InternalResult()

        abcd = self.solve(freqs, **kvargs)
        result.storeSignal('twoport', TwoPort(abcd))

        result.storeSignal('mu', 1/abcd[0,0])
        result.storeSignal('gamma', 1/abcd[0,1])
        result.storeSignal('zeta', 1/abcd[1,0])
        result.storeSignal('beta', 1/abcd[1,1])

        if self.noise:
            inp, inn, outp, outn = self.ports
            
            circuit_voltagesrc = copy(self.c)
            circuit_voltagesrc['VS_TwoPort'] = VS(inp, inn, vac = 1)
            
            circuit_currentsrc = copy(self.c)
            circuit_currentsrc['IS_TwoPort'] = IS(inp, inn, iac = 1)
            
            if self.noise_outquantity == 'i':
                for src in circuit_voltagesrc, circuit_currentsrc:
                    src['VL'] = VS(outp, outn, vac = 0)

            if self.noise_outquantity == 'v':
                res_v = self.NoiseAnalysis(circuit_voltagesrc, inputsrc=circuit_voltagesrc['VS_TwoPort'], \
                                           outputnodes=(outp, outn)).run(freqs)

                res_i = self.NoiseAnalysis(circuit_currentsrc, inputsrc=circuit_currentsrc['IS_TwoPort'], \
                                           outputnodes=(outp, outn)).run(freqs)
            else:
                res_v = self.NoiseAnalysis(circuit_voltagesrc, inputsrc=circuit_voltagesrc['VS_TwoPort'], \
                                           outputsrc=circuit_voltagesrc['VL']).run(freqs)

                res_i = self.NoiseAnalysis(circuit_currentsrc, inputsrc=circuit_currentsrc['IS_TwoPort'], \
                                           outputsrc=circuit_currentsrc['VL']).run(freqs)

            result['Svn'] = res_v['Svninp']
            result['Sin'] = res_i['Sininp']
            
        self.result = result

        return result

    def solve(self, freqs, refnode = gnd, complexfreq = False):
        inp, inn, outp, outn = self.ports
                
        ## Add voltage source at input port and create
        ## copies with output open and shorted respectively
        circuit_vs_open = copy(self.c)

        circuit_vs_open['VS_TwoPort'] = VS(inp, inn, vac=1)

        circuit_vs_shorted = copy(circuit_vs_open)

        circuit_vs_shorted['VL_TwoPort'] = VS(outp, outn, vac=0)

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
