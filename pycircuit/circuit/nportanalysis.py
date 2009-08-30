# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import numpy as np
from nport import *
from pycircuit.circuit import SubCircuit, gnd, R, VS, IS, Branch, circuit
from analysis import Analysis, AC, Noise, TransimpedanceAnalysis, \
    remove_row_col,defaultepar
from pycircuit.post.internalresult import InternalResultDict

np.set_printoptions(precision=4)

class TwoPortAnalysis(Analysis):
    """Analysis to find the 2-ports parameters of a circuit

    The transmission parameters are found as:

    A = v(inp, inn)/v(outp, outn) | io = 0
    B = v(inp, inn)/i(outp, outn) | vo = 0
    C = i(inp, inn)/v(outp, outn) | io = 0
    D = i(inp, inn)/i(outp, outn) | vo = 0

    Examples:

    >>> c = SubCircuit()
    >>> n1 = c.add_node('net1')
    >>> n2 = c.add_node('net2')
    >>> c['R1'] = R(n1, n2, r=9e3)
    >>> c['R2'] = R(n2, gnd, r=1e3)
    >>> res = TwoPortAnalysis(c, n1, gnd, n2, gnd).solve(freqs = np.array([0]))
    >>> print res['mu'].y[0]
    (0.1+0j)
    >>> print res['gamma'].y[0]
    (0.000111111111111+0j)
    >>> print res['zeta'].y[0]
    (1000+0j)
    >>> print res['beta'].y[0]
    (1+0j)

    The transmission parameters are found as:

    A = v(inp, inn)/v(outp, outn) | io = 0
    B = v(inp, inn)/i(outp, outn) | vo = 0
    C = i(inp, inn)/v(outp, outn) | io = 0
    D = i(inp, inn)/i(outp, outn) | vo = 0

    >>> import symbolic; from sympy import simplify, Symbol
    >>> circuit.default_toolkit = symbolic
    >>> c = SubCircuit()
    >>> n1, n2 = c.add_nodes('net1', 'net2')
    >>> c['R1'] = R(n1, n2, r=Symbol('R1',real=True))
    >>> c['R2'] = R(n2, gnd, r=Symbol('R2',real=True))
    >>> symnoise = TwoPortAnalysis(c, n1, gnd, n2, gnd, noise=True, toolkit=symbolic)
    >>> res = symnoise.solve(freqs = np.array([Symbol('s')]), complexfreq=True)
    >>> simplify(res['mu'].y[0])
    R2/(R1 + R2)
    >>> simplify(res['gamma'].y[0])
    1/R1
    >>> simplify(res['zeta'].y[0])
    R2
    >>> simplify(res['beta'].y[0])
    1
    >>> simplify(res['Svn'])
    (4*R1*R2*kT + 4*kT*R1**2)/R2
    >>> simplify(res['Sin'])
    4*kT/R2

    
    """
    
    def __init__(self, circuit, inp, inn, outp, outn, noise = False, 
                 noise_outquantity = 'v', method = 'sparam', 
                 toolkit = None):
        super(TwoPortAnalysis, self).__init__(circuit, toolkit=toolkit)

        self.ports = (inp, inn), (outp, outn)

        self.noise = noise
        self.noise_outquantity = noise_outquantity

        self.method = method
        
    def solve(self, freqs, complexfreq = False, refnode = gnd):
        toolkit = self.toolkit

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
            
            circuit_vs = copy(self.cir)
            circuit_vs['VS_TwoPort'] = VS(inp, inn, vac = 1)
            
            circuit_cs = copy(self.cir)
            circuit_cs['IS_TwoPort'] = IS(inp, inn, iac = 1)
            
            if self.noise_outquantity == 'i':
                for src in circuit_vs, circuit_cs:
                    src['VL'] = VS(outp, outn, vac = 0)

            if self.noise_outquantity == 'v':
                na = Noise(circuit_vs, 
                           inputsrc='VS_TwoPort',
                           outputnodes=(outp, outn), 
                           toolkit=toolkit)
                res_v = na.solve(freqs, complexfreq=complexfreq, 
                                 refnode=refnode)

                na = Noise(circuit_cs, 
                           inputsrc='IS_TwoPort',
                           outputnodes=(outp, outn),
                           toolkit=toolkit)
                res_i = na.solve(freqs, complexfreq=complexfreq,
                                 refnode=refnode)
            else:
                na = Noise(circuit_vs, 
                                inputsrc='VS_TwoPort',
                                outputsrc='VL',
                                toolkit=toolkit)
                res_v = na.solve(freqs, complexfreq=complexfreq,
                                 refnode=refnode)

                na = Noise(circuit_cs, 
                           inputsrc='IS_TwoPort',
                           outputsrc='VL',
                           toolkit=toolkit)
                res_i = na.solve(freqs, complexfreq=complexfreq,
                                 refnode=refnode)

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
        ## b_n = ((2-R0*i_n) - R0*i_n) / 2 / sqrt(R0) = {i_n = (v_n - 2)/R0} =
        ## (2-2*R0*(v_n-2)/R0)/2/sqrt(R0) = (1 - v_n - 2) / sqrt(R0) =
        ## = (v_n - 1) / sqrt(R0)
        ## => S_k_n = b_k / a_n = v_k | k != n
        ## S_n_n = b_n / a_n = v_n - 1
        ##
        ## 
        toolkit = self.toolkit

        # Reference impedance
        r0 = 1

        N = len(self.ports)

        portnumbers = range(N)

        S = np.zeros((N,N), dtype=object)
        
        for n, sourceport in enumerate(self.ports):
            circuit = copy(self.cir)
            
            vs_plus = circuit.add_node('vs_plus')
            
            ## Power source at port n
            circuit['_vs'] = VS(vs_plus, sourceport[1], vac = 2)
            circuit['_rs'] = R(vs_plus, sourceport[0], r = r0)
            
            ## Terminate other ports
            for k, port in enumerate(self.ports):
                if port != sourceport:
                    circuit['_rl%d'%k] = R(port[0], port[1], r = r0)
            
            ## Run AC-analysis
            res = AC(circuit, toolkit=toolkit).solve(freqs, 
                                                     refnode=sourceport[1],
                                                     complexfreq = complexfreq)

            ## Obtain s-parameters
            for k, port in enumerate(self.ports):
                if k == n:
                    S[k,n] = res.v(port[0], port[1]) - 1
                else:
                    S[k,n] = res.v(port[0], port[1])

        ## Find noise wave correlation matrix
        ## The method works as follows:
        ## 1. terminate all ports with z0
        ## 2. calculate transimpedances from a current source in each node
        ##    to the port voltages by using the adjoint Y-matrix.
        ##    As seen above the outgoing wave b is b=V(port_k) / sqrt(z0)
        ##    for a terminated port.
        ## 3. Form a transform matrix T where the columns are the transimpedance
        ##    vectors divided by sqrt(z0)
        ## 4. Calculate the correlation matrix as:
        ##    CS = T * CY * T+
        
        circuit = copy(self.cir)

        ## Terminate ports
        for k, port in enumerate(self.ports):
            circuit['_rl%d'%k] = R(port[0], port[1], r = r0, noisy = False)
        
        ## Calculate transimpedances
        branchlist = [Branch(*port) for port in self.ports]
        
        refnode = self.ports[0][1]
        
        transimpana = TransimpedanceAnalysis(circuit, toolkit = toolkit)

        zmlist = transimpana.solve(freqs,
                                   branchlist,
                                   refnode = self.ports[0][1], 
                                   complexfreq = complexfreq)

        T = np.matrix(zmlist) * r0**-0.5
        
        ## Complex frequency variable
        if complexfreq:
            s = freqs
        else:
            s = 2j*np.pi*freqs

        ## Calculate CY of circuit
        x = np.zeros(circuit.n)
        CY = circuit.CY(x, np.imag(s), epar = self.epar)
        irefnode = circuit.get_node_index(refnode)
        CY, = remove_row_col((CY,), irefnode)

        ## Calculate noise wave correlation matrix
        CS = np.array(T * CY * T.H)

        return NPortS(S, CS, z0=r0)

    def solve_abcd(self, freqs, refnode = gnd, complexfreq = False):
        (inp, inn), (outp, outn) = self.ports
                
        toolkit = self.toolkit

        ## Add voltage source at input port and create
        ## copies with output open and shorted respectively
        circuit_vs_open = copy(self.cir)

        circuit_vs_open['VS_TwoPort'] = VS(inp, inn, vac=1)

        circuit_vs_shorted = copy(circuit_vs_open)

        circuit_vs_shorted['VL_TwoPort'] = VS(outp, outn, vac=0)

        ## Run AC-analysis on the two circuits
        ac_open = AC(circuit_vs_open, toolkit=toolkit)
        ac_shorted = AC(circuit_vs_shorted, toolkit=toolkit)

        res_open = ac_open.solve(freqs, refnode = refnode, 
                                 complexfreq=complexfreq)

        res_shorted = ac_shorted.solve(freqs, refnode = refnode, 
                                     complexfreq=complexfreq)
        
        A = res_open.v(inp, inn) / res_open.v(outp, outn)
        B = res_shorted.v(inp, inn) / res_shorted.i('VL_TwoPort.plus')
        C = res_open.i('VS_TwoPort.minus') / res_open.v(outp, outn)
        D = res_shorted.i('VS_TwoPort.minus') / res_shorted.i('VL_TwoPort.plus')

        return np.array([[A,B],[C,D]], dtype=object)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
