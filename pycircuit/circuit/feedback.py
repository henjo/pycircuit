import numpy as np
from copy import copy
from circuit import SubCircuit, gnd, R, VS, IS, Branch, VCCS, CircuitProxy
from analysis import Analysis, AC, Noise, TransimpedanceAnalysis, \
    remove_row_col,defaultepar,isiterable
from pycircuit.post import InternalResultDict, Waveform

class LoopBreakError(Exception):
    pass

class LoopBreaker(CircuitProxy):
    """Circuit proxy that zeros out dependent sources in the G and C methods"""
    def __init__(self, circuit, inp, inn, outp, outn, parent=None, instance_name=None):
        super(LoopBreaker, self).__init__(circuit, parent, instance_name)
        
        self.inport = tuple([circuit.get_node_index(node) 
                             for node in inp, inn])
        self.outport = tuple([circuit.get_node_index(node) 
                              for node in outp, outn])

    def G(self, x, epar=defaultepar): 
        G = self.device.G(x,epar)

        ## Detect voltage-controlled current sources
        if G[self.outport[0], self.inport[0]] == \
           -G[self.outport[0], self.inport[1]] == \
           -G[self.outport[1], self.inport[0]] == \
           G[self.outport[1], self.inport[1]]:
           for row,col in ((self.outport[0], self.inport[0]),
                           (self.outport[0], self.inport[1]),
                           (self.outport[1], self.inport[0]),
                           (self.outport[1], self.inport[1])):
               G[row,col] = 0
        else:
            raise LoopBreakError('Could not detect dependent source')

        return G

class FeedbackDeviceAnalysis(Analysis):
    """Find loop-gain by replacing a dependent source with a independent one

    Example: 

    >>> cir = SubCircuit()
    >>> cir['M1'] = VCCS('g', 's', gnd, 's', gm = 20e-3)
    >>> cir['RL'] = R('s', gnd)
    >>> cir['VS'] = VS('g', gnd)
 
    >>> ana = FeedbackDeviceAnalysis(cir, 'M1', 'inp', 'inn', 'outp', 'outn')
    >>> res = ana.solve(1e3)
    >>> np.around(res['loopgain'])
    (-20-0j)
    
    """
    def __init__(self, circuit, instance, inp=None, inn=None, outp=None, outn=None, 
                 epar = defaultepar.copy()):
        super(FeedbackDeviceAnalysis, self).__init__(circuit, epar = epar)

        self.inp = inp; self.inn = inn; self.outp = outp; self.outn = outn

        self.device = instance
        
        try:
            circuit[instance]
        except KeyError:
            raise ValueError('Device %s does not exist'%instance)

        
    def solve(self, freqs, complexfreq = False, refnode = gnd):
        circuit_noloop = copy(self.cir)

        circuit_noloop[self.device] = LoopBreaker(circuit_noloop[self.device], 
                                                  self.inp,self.inn,
                                                  self.outp,self.outn,
                                                  parent = circuit_noloop, 
                                                  instance_name = self.device)

        
        x = np.zeros(self.cir.n) ## FIXME, this should come from the DC analysis
        
        epar = self.epar
        G = self.cir.G(x, epar)
        G_noloop = circuit_noloop.G(x, epar)

        C = self.cir.C(x, epar)
        C_noloop = circuit_noloop.C(x, epar)

        ## Refer the voltages to the reference node by removing
        ## the rows and columns that corresponds to this node
        irefnode = self.cir.get_node_index(refnode)
        G,G_noloop,C,C_noloop = remove_row_col((G,G_noloop,C,C_noloop), 
                                               irefnode)
        if complexfreq:
            slist = freqs
        else:
            slist = 2j*np.pi*freqs
        
        ## Calculate return difference
        def Ffunc(s):
            Y = self.toMatrix(G + s*C)
            Y_noloop = self.toMatrix(G_noloop + s*C_noloop)
            return self.det(Y) / self.det(Y_noloop)
        if isiterable(slist):
            F = Waveform(np.array(slist),
                         np.array([Ffunc(s) for s in slist]),
                         xlabels = ('frequency',),
                         xunits = ('Hz',),
                         ylabel = 'F')
        else:
            F = Ffunc(slist)

        ## Calculate loop-gain
        T = F - 1

        result = InternalResultDict()
        result['F'] = F
        result['T'] = T
        result['loopgain'] = -T

        return result
        
if __name__ == "__main__":
    import doctest
    doctest.testmod()
