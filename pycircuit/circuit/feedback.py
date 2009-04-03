import numpy as np
from copy import copy
from pycircuit.circuit import Circuit, SubCircuit, gnd, R, VS, IS, \
    Branch, VCCS, CircuitProxy
from analysis import Analysis, AC, Noise, TransimpedanceAnalysis, \
    remove_row_col,defaultepar,isiterable
from pycircuit.post import InternalResultDict, Waveform
from pycircuit.utilities import combinations

class LoopBreakError(Exception):
    pass

def find_nulling_indices(circuit, x, inp=None, inn=None, outp=None, outn=None, 
                         epar=defaultepar):
    n = circuit.n
    G = circuit.G(x,epar)
    C = circuit.C(x,epar)
    
    terminal_node_indices = [circuit.get_node_index(t) 
                             for t in circuit.terminals]
    if len(circuit.terminals) == 2:
        return [[[terminal_node_indices[0]],[terminal_node_indices[1]]],
                terminal_node_indices]
    else:
        if inp == None:
            allports = combinations(terminal_node_indices, 2)
            all_combinations = combinations(allports, 2)
        else:
            i_inp, i_inn, i_outp, i_outn = [circuit.get_node_index(t)
                                            for t in inp,inn,outp,outn]
            all_combinations = (((i_inp, i_inn), (i_outp, i_outn)),)
        
        found_combinations = set()
        for inport, outport in all_combinations:
            for X in G,C:
                inbranch = Branch(*[circuit.nodes[n] for n in inport])
                outbranch = Branch(*[circuit.nodes[n] for n in outport])

                if outbranch in circuit.branches:
                    ## VCVS
                    i_outbranch = circuit.get_branch_index(outbranch)
                    g = X[i_outbranch, inport[0]]

                    if g == 0:
                        continue
                    
                    stampindex = [i_outbranch, inport]
                    stampfound = np.alltrue(X[stampindex] == np.array([g,-g]))

                    indices = set(range(n)) - set(inport) - set(outport)
                    no_other_crossterms = np.alltrue(X[i_outbranch, 
                                                       list(indices)] == 0)

                    if stampfound and no_other_crossterms:
                        found_combinations.add((inport, outport))
                        null_indices = stampindex
                    

                    ## CCVS
                    if inbranch in circuit.branches:
                        r = X[i_outbranch, i_inbranch]
                    
                        if r == 0:
                            continue

                        indices = set(range(n)) - set([i_inbranch])
                        no_other_crossterms = np.alltrue(X[i_outbranch, 
                                                           list(indices)] == 0)

                        if no_other_crossterms:
                            found_combinations.add((inport, outport))
                            null_indices = [(i_outbranch, i_inbranch)]
                else:
                    ## VCCS
                    gm = X[outport[0], inport[0]]

                    if gm == 0:
                        continue

                    stampindex = [[[outport[0]],[outport[1]]], inport]
                    stampfound = np.alltrue(X[stampindex] == \
                                            np.array([[gm,-gm],[-gm,gm]]))

                    other_cross_rows = []; other_cross_cols = []
                    for row in outport:
                        for col in set(range(n)) - set(inport) - set([row]):
                            other_cross_rows.append(row)
                            other_cross_cols.append(col)

                    no_other_crossterms = np.alltrue(X[other_cross_rows, 
                                                       other_cross_cols] == 0)

                    if stampfound and no_other_crossterms:
                        found_combinations.add((inport, outport))
                        null_indices = stampindex



        if len(found_combinations) == 1:
            return null_indices
            

class LoopBreaker(CircuitProxy):
    """Circuit proxy that zeros out dependent sources in the G and C methods"""
    def __init__(self, circuit, inp, inn, outp, outn, parent=None, 
                 instance_name=None):
        super(LoopBreaker, self).__init__(circuit, parent, instance_name)
        
        self.nulling_indices = find_nulling_indices(circuit, 
                                                    np.zeros(circuit.n),
                                                    inp, inn, outp, outn)
        
        if self.nulling_indices == None:
            raise LoopBreakError('Could not detect dependent source')
        
    def G(self, x, epar=defaultepar): 
        G = self.device.G(x,epar)
        G[self.nulling_indices] = 0
        return G

    def C(self, x, epar=defaultepar): 
        C = self.device.C(x,epar)
        C[self.nulling_indices] = 0
        return C

class FeedbackDeviceAnalysis(Analysis):
    """Find loop-gain by replacing a dependent source with a independent one

    Example: 

    >>> cir = SubCircuit()
    >>> cir['M1'] = VCCS('g', 's', gnd, 's', gm = 20e-3)
    >>> cir['RL'] = R('s', gnd)
    >>> cir['VS'] = VS('g', gnd)
 
    >>> ana = FeedbackDeviceAnalysis(cir, 'M1')
    >>> res = ana.solve(1e3)
    >>> np.around(res['loopgain'])
    (-20-0j)
    
    """
    def __init__(self, circuit, instance, 
                 inp=None, inn=None, outp=None, outn=None, 
                 epar = defaultepar.copy(), toolkit=None):
        super(FeedbackDeviceAnalysis, self).__init__(circuit, epar = epar,
                                                     toolkit=toolkit)

        self.inp = inp; self.inn = inn; self.outp = outp; self.outn = outn

        self.device = instance
        
        try:
            circuit[instance]
        except KeyError:
            raise ValueError('Device %s does not exist'%instance)

        
    def solve(self, freqs, complexfreq = False, refnode = gnd):
        toolkit = self.toolkit

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
            slist = 2j * np.pi * freqs
        
        ## Calculate return difference
        def Ffunc(s):
            Y = toolkit.toMatrix(G + s*C)
            Y_noloop = toolkit.toMatrix(G_noloop + s*C_noloop)
            return toolkit.det(Y) / toolkit.det(Y_noloop)
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


class LoopVS(VS):
    def u(self, t=0.0, epar=defaultepar, analysis=None):
        if analysis == 'feedback':
            return np.array([0, 0, -self.ipar.vac], dtype=object)
        else:
            return np.zeros(self.n)

class LoopIS(IS):
    def u(self, t=0.0, epar=defaultepar, analysis=None):
        if analysis == 'feedback':
            return np.array([self.ipar.iac, -self.ipar.iac])
        else:
            return np.zeros(self.n)

class LoopProbe(SubCircuit):
    terminals = ['inp', 'inn', 'outp', 'outn']

    def __init__(self, *args, **kvargs):
        super(LoopProbe, self).__init__(*args, **kvargs)

        self['vinj'] = LoopVS(self.nodenames['inp'], self.nodenames['outp'])
        self['iinj'] = LoopIS(self.nodenames['inp'], self.nodenames['inn'])

def find_loopprobe(circuit, path=[]):
    for name, e in circuit.elements.items():
        if isinstance(e, LoopProbe):
            yield ('.'.join(path + [name]), e)
        if isinstance(e, SubCircuit):
            find_loopprobe(e, path + [name])
                
        
class FeedbackLoopAnalysis(Analysis):
    """Find loop-gain by breaking all loops with a LoopProbe element
    """
    def __init__(self, circuit, epar = defaultepar.copy(), toolkit=None):
        super(FeedbackLoopAnalysis, self).__init__(circuit, epar = epar, 
                                                   toolkit=toolkit)
        
        ## Find LoopProbe instance
        loopprobes = list(find_loopprobe(circuit))
        
        if len(loopprobes) != 1:
            raise ValueError('The circuit must contain exactly one LoopProbe instance')
        
        self.loopprobe_name, self.loopprobe = loopprobes[0]


    def solve(self, freqs, refnode = gnd, complexfreq = False):
        toolkit = self.toolkit

        x = np.zeros(self.cir.n) ## FIXME, this should be replaced by DC-analysis
        self.loopprobe['vinj'].ipar.vac = 1 
        self.loopprobe['iinj'].ipar.iac = 0 

        ac_vinj = AC(self.cir, toolkit=toolkit)

        res_vinj = ac_vinj.solve(freqs, refnode = refnode, 
                                 complexfreq=complexfreq,
                                 u = self.cir.u(x, analysis='feedback'))

        self.loopprobe['vinj'].ipar.vac = 0 
        self.loopprobe['iinj'].ipar.iac = 1 

        ac_iinj = AC(self.cir, toolkit=toolkit)

        res_iinj = ac_iinj.solve(freqs, refnode = refnode, 
                                 complexfreq=complexfreq,
                                 u = self.cir.u(x, analysis='feedback'))

        B = res_vinj.i(self.loopprobe_name + '.vinj.plus')
        D = res_vinj.v(self.loopprobe_name + '.inp', self.loopprobe_name + '.inn')
        A = res_iinj.i(self.loopprobe_name + '.vinj.plus')
        C = res_iinj.v(self.loopprobe_name + '.inp', self.loopprobe_name + '.inn')

        T = (2*(A*D - B*C) - A + D) / (2*(B*C - A*D) + A - D + 1)
        
        result = InternalResultDict()
        result['F'] = 1 + T
        result['T'] = T
        result['loopgain'] = -T
        
        return result

if __name__ == "__main__":
    import doctest
    doctest.testmod()
