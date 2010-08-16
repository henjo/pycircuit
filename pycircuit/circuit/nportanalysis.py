# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

import numpy as np
from nport import *
from pycircuit.utilities import isiterable
from pycircuit.circuit import SubCircuit, gnd, R, VS, IS, Branch, circuit
from analysis import Analysis, remove_row_col,defaultepar
from analysis_ss import AC, Noise, TransimpedanceAnalysis, dc_steady_state

from pycircuit.post.internalresult import InternalResultDict

np.set_printoptions(precision=4)

class ISInternal(IS):
    """Current source that only responds to 'internalac' analysis"""
    def u(self, t=0, epar=defaultepar, analysis=None):
        if analysis == 'internalac':
            return self.toolkit.array([self.ipar.iac, -self.ipar.iac])
        else:
            return self.toolkit.array([0, 0])

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
                           toolkit=toolkit,
                           epar=self.epar)
                res_v = na.solve(freqs, complexfreq=complexfreq, 
                                 refnode=refnode)

                na = Noise(circuit_cs, 
                           inputsrc='IS_TwoPort',
                           outputnodes=(outp, outn),
                           toolkit=toolkit,
                           epar=self.epar)
                res_i = na.solve(freqs, complexfreq=complexfreq,
                                 refnode=refnode)
            else:
                na = Noise(circuit_vs, 
                           inputsrc='VS_TwoPort',
                           outputsrc='VL',
                           toolkit=toolkit,
                           epar=self.epar)
                res_v = na.solve(freqs, complexfreq=complexfreq,
                                 refnode=refnode)

                na = Noise(circuit_cs, 
                           inputsrc='IS_TwoPort',
                           outputsrc='VL',
                           toolkit=toolkit,
                           epar=self.epar)
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
        
        circuit = copy(self.cir)

        refnode = self.ports[0][1]

        ## Place power sources at the ports
        for n, sourceport in enumerate(self.ports):
            circuit['_is%d'%n] = ISInternal(sourceport[1], sourceport[0], iac = 0, toolkit=toolkit)
            circuit['_rl%d'%n] = R(sourceport[1], sourceport[0], r = r0, noisy=False, toolkit=toolkit)

        if toolkit.symbolic:
            ## For now just one frequency is allowed
            assert not isiterable(freqs)
            
            G, C, u, x, ss = dc_steady_state(circuit, freqs, 
                                             refnode,
                                             toolkit,
                                             complexfreq = complexfreq,
                                             epar = self.epar)
            Y = C * ss + G
            detY =  toolkit.det(Y)

        for n, sourceport in enumerate(self.ports):
            ## Add stimulus to the port
            circuit['_is%d'%n].ipar.iac = 2 / r0

            ## If symbolic the s-parameters are calculated using co-factors
            if toolkit.symbolic:
                u_wrefnode = circuit.u(x, analysis='internalac', epar=self.epar)
                (u,) = circuit.remove_refnode((u_wrefnode,), refnode)
                ## Calculate s-parameters using cofactors
                for k, port in enumerate(self.ports):
                    resname = "v(%s,%s)"%(port[0], port[1])

                    res = linearsolver_partial(Y, u, refnode, [resname],
                                               circuit, toolkit, detY=detY)
                    if k == n:
                        S[k,n] = res[resname] - 1
                    else:
                        S[k,n] = res[resname]

            else:
                ## Run AC-analysis
                acana = AC(circuit, toolkit=toolkit, analysis='internalac')
                res = acana.solve(freqs, refnode=refnode,
                                  complexfreq = complexfreq)
                ## Obtain s-parameters
                for k, port in enumerate(self.ports):
                    if k == n:
                        S[k,n] = res.v(port[0], port[1]) - 1
                    else:
                        S[k,n] = res.v(port[0], port[1])

            ## Clear stimulus to the port
            circuit['_is%d'%n].ipar.iac = 0

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

        ## Calculate transimpedances
        branchlist = [Branch(*port) for port in self.ports]
        
        transimpana = TransimpedanceAnalysis(circuit, toolkit=toolkit)

        zmlist = transimpana.solve(freqs, branchlist, refnode=refnode,
                                   complexfreq=complexfreq)

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

def linearsolver_partial(Y, u, refnode, selected_res, cir, toolkit, detY=None):
    """Solve linear system Y * x + u = 0 and return a dictionary of selected result

    The function should only be used for symbolic calculations since more
    efficient methods exists for numeric problems.

    The selected_res argument is a list/tuple of desired results in the form 
    "v(nodea, nodeb)" or "v(nodea)" for voltage potentials.

    """
    if detY == None:
        detY = toolkit.det(Y)
    
    uindices = toolkit.nonzero(u)

    result = {}

    for res_str in selected_res:
        if res_str[0:2] != 'v(' or res_str[-1] != ')':
            raise ValueError('Invalid result selector: %s'%res_str)

        nodes = res_str[2:-1].split(',')

        nodes_indices = [cir.get_node_index(node, refnode) for node in nodes]

        ## xj = (sum_i wi * cofactor(Y,i,j)) / det Y
        if Y.shape == (1,1):
            num = -1
        else:
            num = 0
            for ui in uindices:
                for sign, nodeindex in zip([1,-1], nodes_indices):
                    if nodeindex != None:
                        num += sign * -u[ui] * toolkit.cofactor(Y, ui, nodeindex)
                        
        result[res_str] = num / detY

    return result

if __name__ == "__main__":
    import doctest
    doctest.testmod()
