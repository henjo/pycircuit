# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import numpy as npy
from circuit import SubCircuit, gnd, R, VS, IS
from analysis import Analysis, AC, Noise
from pycircuit.post.internalresult import InternalResultDict
from copy import copy
import constants

class NPort(object):
    """Class that represents an n-port with optional noise parameters

    Attributes
    ----------
    n -- number of ports
    passive -- True if n-port is passive
    noise -- True if n-port has noise parameters 
    S -- S-parameter matrix
    Y -- Y-parameter matrix
    Z -- Z-parameter matrix
    A -- ABCD-parameter matrix
    CS -- Noise wave correlation matrix
    CY -- Y-parameter noise correlation matrix
    CZ -- Z-parameter noise correlation matrix
    CA -- ABCD-parameter noise correlation matrix

    """
    passive = False

    def __mul__(self, a):
        """Cascade of two n-ports"""
        anport = NPortA(dot(self.A, a.A), dot(dot(self.A, a.CA), npy.conjugate(self.A).T) + self.CA )
        return self.__class__(anport)

    def __floordiv__(self, a):
        """Parallel of two n-ports

        >>> import sympy as S
        >>> a,b,c,d = S.symbols('abcd')
        >>> A = NPortA(npy.array([[a,b], [c,d]]))
        >>> (A // A).A
        array([[a, 0.5*b],
           [-2c, d]], dtype=object)

        """
        ynport = NPortY(self.Y + a.Y, self.CY + a.CY)
        return self.__class__(ynport)

    def series(self, a):
        """Series connection with another n-port

        >>> import sympy as S
        >>> a,b,c,d = S.symbols('abcd')
        >>> A = NPortA(npy.array([[a,b], [c,d]]))
        >>> A.series(A).A
        [[ a 2*b]
         [ c/2 d]]
        
        """
        znport = NPortZ(self.Z + a.Z, self.CZ + a.CZ)
        
        return self.__class__(znport)

    def noisy_passive_nport(self, T=290):
        """Returns an n-port with noise parameters set if passive"""

        if not self.passive:
            raise ValueError("Cannot calculate noise-correlation matrix of non-passive n-port")
            
        ynport = NPortY(self)
        
        ynport.CY = 2 * constants.kboltzmann * T * npy.real(ynport.Y)

        return ynport        

class NPortY(NPort):
    """Two-port class where the internal representation is the Y-parameters"""
    
    def __init__(self, Y, CY = None):
        if isinstance(Y, NPort):
            self.Y = Y.Y
        else:
            self.Y = Y
        
        if CY == None:
            self.CY = npy.zeros(npy.shape(self.Y))
        else:
            self.CY = CY

        self.n = npy.size(self.Y,0)
        
    @property
    def A(self):
        """Return chain parameters (ABCD)"""
        if self.n != 2:
            raise ValueError('N-port must be a 2-port')

        Y = self.Y
        d = Y[0,0] * Y[1,1] - Y[0,1] * Y[1,0]
        return npy.array([[-Y[1,1] / Y[1,0], -1.0 / Y[1,0]],
                          [-d / Y[1,0], -Y[0,0] / Y[1,0]]])

    @property
    def Z(self):
        """Return Z-parameter matrix"""
        return npy.invert(self.Y)

    @property
    def S(self, z0 = 50.0):
        """Return scattering parameters"""
        
        E = npy.mat(npy.eye(self.n, self.n))
        Zref = z0 * E
        Gref = 1 / npy.sqrt(npy.real(z0)) * E
        return Gref * (E - Zref * Y) * (E + Zref * Y)**-1 * Gref**-1

    @property
    def CZ(self):
        Z = mat(self.Z)
        return npy.array(Z * self.CY * Z.H)

    @property
    def CS(self):
        E = npy.mat(npy.eye(self.n, self.n))
        return npy.array((E + S) * sympy.mat(self.CY) * (E + S).H / 4)

    @property
    def CA(self):
        T = npy.mat(self.A)
        T[0,0] = 0
        T[1,0] = 1

        return npy.array(T * npy.mat(self.CY) * T.H)

class NPortZ(NPort):
    """Two-port class where the internal representation is the Z-parameters"""
    
    def __init__(self, Z, CZ = None):
        if isinstance(Z, NPort):
            self.Z = Z.Z
        else:
            self.Z = Z
        
        if CZ == None:
            self.CZ = npy.zeros(npy.shape(self.Y))
        else:
            self.CZ = CZ

        self.n = npy.size(self.Z,0)

    @property
    def A(self):
        """Return chain parameters (ABCD)"""
        if self.n != 2:
            raise ValueError('N-port must be a 2-port')

        Z = self.Z
        d = Z[0,0] * Z[1,1] - Z[0,1] * Z[1,0]
        return npy.array([[-Z[0,0] / Z[1,0], -d / Z[1,0]],
                          [1.0 / Z[1,0], Z[1,1] / Z[1,0]]])
    
    @property
    def Y(self):
        """Return Z-parameter matrix"""
        return npy.invert(self.Z)

    @property
    def S(self, z0 = 50.0):
        """Return scattering parameters"""
        Z = self.Z
        E = npy.mat(npy.eye(self.n, self.n))
        Zref = z0 * E
        Gref = 1 / npy.sqrt(npy.real(z0)) * E
        return Gref * (Z - Zref) * (Z + Zref)**-1 * Gref**-1

    @property
    def CY(self):
        Y = mat(self.Y)
        return Y * self.CZ * Y.H

    @property
    def CS(self):
        E = npy.mat(npy.eye(self.n, self.n))
        return (E - S) * npy.mat(self.CZ) * (E - S).H / 4

    @property
    def CA(self):
        T = npy.matrix([[1, -self.A[0,0]], [0, -self.A[1,0]]])

        return T * npy.mat(self.CZ) * T.H

class NPortA(NPort):
    """Two-port class where the internal representation is the ABCD-parameters"""

    def __init__(self, A, CA = None):
        if isinstance(A, NPort):
            self.A = A.A
        else:
            self.A = A

        if npy.shape(self.A) != (2,2):
            raise ValueError('Can only create ABCD-two ports')
        
        if CA == None:
            self.CA = npy.zeros(npy.shape(self.Y))
        else:
            self.CA = CA

        self.n = 2

    @property
    def Z(self):
        """Return Z-parameter matrix"""
        A = self.A
        d = A[0,0] * A[1,1] - A[0,1] * A[1,0]

        return npy.array([[A[0,0] / A[1,0], d / A[1,0]],
                        [1.0 / A[1,0], A[1,1] / A[1,0]]])
    

    @property
    def Y(self):
        """Return Y-parameter matrix"""
        A = self.A
        d = A[0,0] * A[1,1] - A[0,1]*A[1,0]

        return npy.array([[A[1,1] / A[0,1], -d / A[0,1]],
                        [-1.0 / A[0,1], A[0,0] / A[0,1]]])
    
    @property
    def S(self, z0 = 50.0):
        """Return scattering parameters
        
        >>> abcd = npy.array([[  5.90000000e-01,   8.05000000e+01], \
                              [  4.20000000e-03,   1.59000000e+00]])
        >>> P = NPortA(abcd)
        >>> P.S
        array([[ 0.1,  0.3,  0.5,  0.6]])

        >>> 
        """
        a,b,c,d = self.A[0,0], self.A[0,1], self.A[1,0], self.A[1,1]

        A = npy.array([[a + b / z0 - c * z0 - d, 2 * (a * d - b * c),
                        2,                       -a+b/z0-c*z0+d]])
        return 1/(a + b / z0 + c * z0 + d) * A

    @property
    def CY(self):
        Y = self.Y
        T = npy.matrix([[-Y[0,0], 1], [-Y[1,0], 0]])

        return npy.array(T * npy.mat(self.CA) * T.H)

    @property
    def CZ(self):
        Z = self.Z
        T = npy.matrix([[1, -Z[0,0]], [0, -Z[1,0]]])

        return npy.array(T * npy.mat(self.CA) * T.H)

    @property
    def CS(self):
        E = npy.mat(npy.eye(self.n, self.n))
        return npy.array((E - S) * npy.mat(self.CZ) * (E - S).H / 4)

    def __str__(self):
        return self.__class__.__name__ + '(' + repr(self.A) + ')'


class NPortS(NPort):
    """Two-port class where the internal representation is the S-parameters"""
    
    def __init__(self, S, CS = None, z0 = 50):
        self.z0 = z0
        
        if isinstance(S, NPort):
            self.S = S.S
        else:
            self.S = S
        
        if CS == None:
            self.CS = npy.zeros(npy.shape(self.S))
        else:
            self.CS = CS

        self.n = npy.size(self.S,0)

    @property
    def A(self):
        """Return chain parameters (ABCD)
        
        >>> S = npy.array([[0.1,0.3],[0.5,0.6]])
        >>> NPortS(S).A
        array([[0.59, 80.5],
               [0.0, 1.59]], dtype=object)

        """
        s = self.S
        z0 = self.z0
        
        a = ((1 + s[0,0]) * (1 - s[1,1]) + s[0,1]*s[1,0]) / (2 * s[1,0])
        b = z0 * ((1 + s[0,0]) * (1 + s[1,1]) - s[0,1]*s[1,0]) / (2 * s[1,0])
        c = 1 / z0 * ((1 - s[0,0]) * (1 - s[1,1]) - s[0,1]*s[1,0]) / (2 * s[1,0])
        d = ((1 - s[0,0]) * (1 + s[1,1]) + s[0,1]*s[1,0]) / (2 * s[1,0])
        
        return npy.array([[a,b],[c,d]], dtype=object)

    @property
    def Z(self):
        """Return Z-parameter matrix"""
        S = npy.mat(self.S)
        E = npy.mat(npy.eye(self.n, self.n))
        Zref = self.z0 * E
        Gref = 1 / npy.sqrt(npy.real(self.z0)) * E
        return npy.array(Gref**-1 * (E - S)**-1 * (S + E) * Zref * Gref)

    @property
    def Y(self):
        """Return Z-parameter matrix"""
        S = npy.mat(self.S)
        E = npy.mat(npy.eye(self.n, self.n))
        Zref = self.z0 * E
        Gref = 1 / npy.sqrt(npy.real(self.z0)) * E
        return npy.array(Gref**-1 * Zref**-1 * (S + E)**-1 * (E - S) * Gref)

    @property
    def CY(self):
        Y = self.Y
        E = npy.mat(npy.eye(self.n, self.n))
        T = E + Y
        return T * npy.mat(self.CS) * T.H

    @property
    def CZ(self):
        Z = self.Z
        E = npy.mat(npy.eye(self.n, self.n))
        T = E + Z
        return T * npy.mat(self.CS) * T.H

    @property
    def CA(self):
        T = npy.matrix([[1, -self.A[0,0]], [0, -self.A[1,0]]])

        return T * npy.mat(self.CZ) * T.H

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
            abcd = result['twoport'].A
        else:
            abcd = self.solve_abcd(freqs, refnode=refnode, complexfreq=complexfreq)
            result['twoport'] = NPortA(abcd)

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
        >>> mu = 1/twoport.A[0,0]
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
        ##
        ## 
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

        return NPortS(S, z0=r0)

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
