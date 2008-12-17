# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import numpy as npy
from circuit import SubCircuit, gnd, R, VS, IS
from analysis import Analysis, AC, Noise
from pycircuit.post.internalresult import InternalResultDict
from copy import copy

class NPort(object):
    """Class to perform calculations and conversion of n-port circuits

    Attributes
    ----------
    n -- number of ports
    K -- (numpy array) transmission parameter matrix, the size is n x n

    """

    def __init__(self, K):
        assert npy.size(K,0) == npy.size(K,1)
        
        self.K = K

        self.n = npy.size(K,0)

    @classmethod
    def from_abcd(cls, abcd):
        """Create an NPort from ABCD parameters"""
        return cls(abcd)
        
    @classmethod
    def from_y(cls, Y):
        """Create an NPort from Y-parameters"""
        d = Y[0,0] * Y[1,1] - Y[0,1] * Y[1,0]
        return cls(npy.array([[-Y[1,1] / Y[1,0], -1.0 / Y[1,0]],
                              [-d / Y[1,0], -Y[0,0] / Y[1,0]]]))

    @classmethod
    def from_z(cls, Z):
        """Create an NPort from Z-parameters"""
        d = Z[0,0] * Z[1,1] - Z[0,1] * Z[1,0]
        return cls(npy.array([[-Z[0,0] / Z[1,0], -d / Z[1,0]],
                              [1.0 / Z[1,0], Z[1,1] / Z[1,0]]]))

    @classmethod
    def from_s(cls, s, z0=50.0):
        """Create an NPort from scattering parameters
        
        >>> S = npy.array([[0.1,0.3],[0.5,0.6]])
        >>> print TwoPort.from_s(S)
        TwoPort(array([[  5.90000000e-01,   8.05000000e+01],
               [  4.20000000e-03,   1.59000000e+00]]))
        """

        a = ((1 + s[0,0]) * (1 - s[1,1]) + s[0,1]*s[1,0]) / (2 * s[1,0])
        b = z0 * ((1 + s[0,0]) * (1 + s[1,1]) - s[0,1]*s[1,0]) / (2 * s[1,0])
        c = 1 / z0 * ((1 - s[0,0]) * (1 - s[1,1]) - s[0,1]*s[1,0]) / (2 * s[1,0])
        d = ((1 - s[0,0]) * (1 + s[1,1]) + s[0,1]*s[1,0]) / (2 * s[1,0])
        
        return cls(npy.array([[a,b],[c,d]], dtype=object))

    def __mul__(self, a):
        """Cascade of two n-ports"""
        return NPort(dot(self.K, a.k))

    def __floordiv__(self, a):
        """Parallel of two n-ports

        >>> import sympy as S
        >>> a,b,c,d = S.symbols('abcd')
        >>> A = TwoPort(npy.array([[a,b], [c,d]]))
        >>> print A // A
        array([[a, 0.5*b],
           [-2c, d]], dtype=object)

        """
        return NPort.from_y(self.Y + a.Y)

    def series(self, a):
        """Series connection with another n-port

        >>> import sympy as S
        >>> a,b,c,d = S.symbols('abcd')
        >>> A = TwoPort(npy.array([[a,b], [c,d]]))
        >>> A.series(A).K
        [[ a 2*b]
         [ c/2 d]]
        
        """
        return NPort.from_z(self.Z + a.Z)

    @property
    def abcd(self):
        """Return chain parameters (ABCD)"""
        return self.K

    @property
    def Z(self):
        """Return Z-parameter matrix"""
        if self.n != 2:
            raise Exception('Cannot handle %d-port parameters yet'%self.n)
        
        A = self.K
        d = A[0,0] * A[1,1] - A[0,1] * A[1,0]

        return npy.array([[A[0,0] / A[1,0], d / A[1,0]],
                        [1.0 / A[1,0], A[1,1] / A[1,0]]])
    

    @property
    def Y(self):
        """Return Y-parameter matrix"""
        if self.n != 2:
            raise Exception('Cannot handle %d-port parameters yet'%self.n)

        A = self.K
        d = A[0,0] * A[1,1] - A[0,1]*A[1,0]

        return npy.array([[A[1,1] / A[0,1], -d / A[0,1]],
                        [-1.0 / A[0,1], A[0,0] / A[0,1]]])
    

    @property
    def S(self, z0 = 50.0):
        """Return scattering parameters
        
        >>> abcd = npy.array([[  5.90000000e-01,   8.05000000e+01], \
                              [  4.20000000e-03,   1.59000000e+00]])
        >>> P = NPort.from_abcd(abcd)
        >>> P.S
        array([[ 0.1,  0.3,  0.5,  0.6]])

        >>> 
        """
        if self.n != 2:
            raise Exception('Cannot handle %d-port parameters yet'%self.n)

        a,b,c,d = self.K[0,0], self.K[0,1], self.K[1,0], self.K[1,1]

        A = npy.array([[a + b / z0 - c * z0 - d, 2 * (a * d - b * c),
                        2,                       -a+b/z0-c*z0+d]])
        return 1/(a + b / z0 + c * z0 + d) * A

    def __str__(self):
        return self.__class__.__name__ + '(' + repr(self.K) + ')'

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
    >>> n1 = c.add_node('net1')
    >>> n2 = c.add_node('net2')
    >>> c['R1'] = R(n1, n2, r=9e3)
    >>> c['R2'] = R(n2, gnd, r=1e3)
    >>> res = TwoPortAnalysis(c, n1, gnd, n2, gnd).solve(freqs = npy.array([0]))
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
    
    def __init__(self, circuit, inp, inn, outp, outn, noise = False, 
                 noise_outquantity = 'v', method = 'sparam'):
        self.c = circuit

        self.ports = (inp, inn), (outp, outn)

        self.noise = noise
        self.noise_outquantity = noise_outquantity

        self.method = 'sparam'
        
    def solve(self, freqs, complexfreq = False, refnode = gnd):
        result = InternalResultDict()

        if self.method == 'sparam':
            result['twoport'] = self.solve_s(freqs, complexfreq=complexfreq)
            abcd = result['twoport'].abcd
        else:
            abcd = self.solve_abcd(freqs, refnode=refnode, complexfreq=complexfreq)
            result['twoport'] = TwoPort(abcd)

        result['mu'] = 1 / abcd[0,0]
        result['gamma'] = 1 / abcd[0,1]
        result['zeta'] = 1 / abcd[1,0]
        result['beta'] = 1 / abcd[1,1]

        if self.noise:
            (inp, inn), (outp, outn) = self.ports
            
            circuit_vs = copy(self.c)
            circuit_vs['VS_TwoPort'] = VS(inp, inn, vac = 1)
            
            circuit_cs = copy(self.c)
            circuit_cs['IS_TwoPort'] = IS(inp, inn, iac = 1)
            
            if self.noise_outquantity == 'i':
                for src in circuit_vs, circuit_cs:
                    src['VL'] = VS(outp, outn, vac = 0)

            if self.noise_outquantity == 'v':
                res_v = self.NoiseAnalysis(circuit_vs, 
                                           inputsrc=circuit_vs['VS_TwoPort'],
                                           outputnodes=(outp, outn)
                                           ).solve(freqs, complexfreq=complexfreq)

                res_i = self.NoiseAnalysis(circuit_cs, 
                                           inputsrc=circuit_cs['IS_TwoPort'],
                                           outputnodes=(outp, outn)
                                           ).solve(freqs, complexfreq=complexfreq)
            else:
                res_v = self.NoiseAnalysis(circuit_vs, 
                                           inputsrc=circuit_vs['VS_TwoPort'],
                                           outputsrc=circuit_vs['VL']
                                           ).solve(freqs, complexfreq=complexfreq)

                res_i = self.NoiseAnalysis(circuit_cs, 
                                           inputsrc=circuit_cs['IS_TwoPort'],
                                           outputsrc=circuit_cs['VL']
                                           ).solve(freqs, complexfreq=complexfreq)

            result['Svn'] = res_v['Svninp']
            result['Sin'] = res_i['Sininp']
            
        self.result = result

        return result

    def solve_s(self, freqs, complexfreq = False):
        """Calculate scattering (s) parameters of circuit

        >>> c = SubCircuit()
        >>> n1 = c.add_node('net1')
        >>> n2 = c.add_node('net2')
        >>> c['R1'] = R(n1, n2, r=9e3)
        >>> c['R2'] = R(n2, gnd, r=1e3)


        >>> an = TwoPortAnalysis(c, n1, gnd, n2, gnd)
        >>> twoport = an.solve_s(freqs = 0)
        >>> mu = 1/twoport.abcd[0,0]
        >>> print mu
        (0.1+0j)

        """

        ## The s-parameters of an n-port are defined as
        ## B = S * A
        ## where A and B are column-vectors of size n of 
        ## the ingoing and outgoing waves respectively
        ## S is an NxN matrix containing the s-parameters
        ## The elements of the A and B vectors are defined as:
        ## 
        ## a_n = (v_n + Z0 * i_n) / 2 * 1 / sqrt(|Re Z0|)
        ## b_n = (v_n - Z0 * i_n) / 2 * 1 / sqrt(|Re Z0|)
        ## 
        ## where v is the port-voltage and i the current flowing 
        ## into the device, Z0 is an arbitrary chosen impedance
        ##
        ## The straight-forward method to calculate the S matrix
        ## that is used here is to run n AC-analyses and in each
        ## analysis connect a voltage source of voltage 2V
        ## in series with a resistor with resistance R0 ohm
        ## (The chosen Z0 value is real). The other ports
        ## are terminated with resistance R0.
        ## The s-parameters S_k_n can now be calculated as:
        ## S_k_n = b_k / a_n
        ## where
        ## a_n = ((2-R0*i_n) + R0 * i+n) / 2 / sqrt(R0) = 1 / sqrt(R0)
        ## b_k = (v_k + v_k) / 2 / sqrt(R0) = v_k / sqrt(R0) | k != n
        ## b_n = ((2-R0*i_n) - R0*i_n) / 2 / sqrt(R0) = {i_n = (v_n - 2)/R0} = (2-2*R0*(v_n-2)/R0)/2/sqrt(R0) =
        ##     = (1 - v_n - 2) / sqrt(R0) = (v_n - 1) / sqrt(R0)
        ## => S_k_n = b_k / a_n = v_k | k != n
        ## S_n_n = b_n / a_n = v_n - 1
        
        # Reference impedance
        import sympy
        r0 = 1

        N = len(self.ports)

        portnumbers = range(N)

        S = npy.zeros((N,N), dtype=object)
        
        for n, sourceport in enumerate(self.ports):
            circuit = copy(self.c)
            
            vs_plus = circuit.add_node('vs_plus')
            
            ## Power source at port n
            circuit['_vs'] = VS(vs_plus, sourceport[1], vac = 2)
            circuit['_rs'] = R(vs_plus, sourceport[0], r = r0)
            
            ## Terminate other ports
            for k, port in enumerate(self.ports):
                if port != sourceport:
                    circuit['_rl%d'%k] = R(port[0], port[1], r = r0)
            
            ## Run AC-analysis
            res = self.ACAnalysis(circuit).solve(freqs, refnode=sourceport[1],
                                                 complexfreq = complexfreq)

            ## Obtain s-parameters
            for k, port in enumerate(self.ports):
                if k == n:
                    S[k,n] = res.v(port[0], port[1]) - 1
                else:
                    S[k,n] = res.v(port[0], port[1])

        return TwoPort.from_s(S, z0=r0)

    def solve_abcd(self, freqs, refnode = gnd, complexfreq = False):
        (inp, inn), (outp, outn) = self.ports
                
        ## Add voltage source at input port and create
        ## copies with output open and shorted respectively
        circuit_vs_open = copy(self.c)

        circuit_vs_open['VS_TwoPort'] = VS(inp, inn, vac=1)

        circuit_vs_shorted = copy(circuit_vs_open)

        circuit_vs_shorted['VL_TwoPort'] = VS(outp, outn, vac=0)

        ## Run AC-analysis on the two circuits
        ac_open = self.ACAnalysis(circuit_vs_open)
        ac_shorted = self.ACAnalysis(circuit_vs_shorted)

        res_open = ac_open.solve(freqs, refnode = refnode, 
                                 complexfreq=complexfreq)

        res_shorted = ac_shorted.solve(freqs, refnode = refnode, 
                                     complexfreq=complexfreq)
        
        A = res_open.v(inp, inn) / res_open.v(outp, outn)
        B = res_shorted.v(inp, inn) / res_shorted.i('VL_TwoPort.plus')
        C = res_open.i('VS_TwoPort.minus') / res_open.v(outp, outn)
        D = res_shorted.i('VS_TwoPort.minus') / res_shorted.i('VL_TwoPort.plus')

        return npy.array([[A,B],[C,D]], dtype=object)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
