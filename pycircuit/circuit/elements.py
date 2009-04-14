# -*- coding: latin-1 -*-
# Copyright (c) 2008 Pycircuit Development Team
# See LICENSE for details.

from circuit import *

class R(Circuit):
    """Resistor element

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> c['R']
    R('plus','minus',r=1000.0,noisy=True)
    >>> c.G(zeros(2))
    array([[0.001, -0.001],
           [-0.001, 0.001]], dtype=object)
    >>> c = SubCircuit()
    >>> n2=c.add_node('2')
    >>> c['R'] = R(n1, n2, r=1e3)
    >>> c.G(zeros(2))
    array([[0.001, -0.001],
           [-0.001, 0.001]], dtype=object)

    """
    terminals = ['plus', 'minus']
    instparams = [Parameter(name='r', desc='Resistance', unit='ohm', 
                            default=1e3),
                  Parameter(name='noisy', desc='No noise', unit='', 
                            default=True),
                  ]

    def G(self, x, epar=defaultepar):
        g = 1/self.ipar.r
        return  array([[g, -g],
                        [-g, g]], dtype=object)

    def CY(self, x, w, epar=defaultepar):
        if self.ipar.noisy:
            if 'kT' in epar:
                iPSD = 4*epar.kT / self.ipar.r
            else:
                iPSD = 4*kboltzmann * epar.T / self.ipar.r
        else:
            iPSD = 0

        return  array([[iPSD, -iPSD],
                       [-iPSD, iPSD]], dtype=object)

class G(Circuit):
    """Conductor element

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['G'] = G(n1, gnd, g=1e-3)
    >>> c['G']
    G('plus','minus',g=0.001,nonoise=False)
    >>> c.G(zeros(2))
    array([[0.001, -0.001],
           [-0.001, 0.001]], dtype=object)

    """
    terminals = ['plus', 'minus']
    instparams = [Parameter(name='g', desc='Conductance', unit='S', 
                            default=1e-3),
                  Parameter(name='nonoise', 
                            desc='If true the conductance is noiseless', 
                            unit='', default=False),
                  ]

    def G(self, x, epar=defaultepar):
        g = self.ipar.g
        return  array([[g, -g],
                        [-g, g]], dtype=object)

    def CY(self, x, w, epar=defaultepar):
        if not self.ipar.nonoise:
            if 'kT' in epar:
                iPSD = 4*epar.kT*self.ipar.g
            else:
                iPSD = 4*kboltzmann*epar.T*self.ipar.g
            return  array([[iPSD, -iPSD],
                           [-iPSD, iPSD]], dtype=object)
        else:
            return super(G, self).CY(x, w, epar=epar)
        

class C(Circuit):
    """Capacitor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['C'] = C(n1, gnd, c=1e-12)
    >>> c.G(zeros(2))
    array([[0, 0],
           [0, 0]], dtype=object)
    >>> c.C(zeros(2))
    array([[1e-12, -1e-12],
           [-1e-12, 1e-12]], dtype=object)

    """

    terminals = ['plus', 'minus']    
    instparams = [Parameter(name='c', desc='Capacitance', 
                            unit='F', default=1e-12)]

    def C(self, x, epar=defaultepar):
        return array([[self.ipar.c, -self.ipar.c],
                      [-self.ipar.c, self.ipar.c]])

class L(Circuit):
    """Inductor

    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['L'] = L(n1, gnd, L=1e-9)
    >>> c.G(zeros(3))
    array([[0.0, 0.0, 1.0],
           [0.0, 0.0, -1.0],
           [1.0, -1.0, 0.0]], dtype=object)
    >>> c.C(zeros(3))
    array([[0, 0, 0],
           [0, 0, 0],
           [0, 0, 1e-09]], dtype=object)
    """
    terminals = ['plus', 'minus']    
    instparams = [Parameter(name='L', desc='Inductance', 
                            unit='H', default=1e-9)]

    _G = array([[0.0 , 0.0, 1.0],
                [0.0 , 0.0, -1.0],
                [1.0 , -1.0, 0.0]])
    def __init__(self, plus, minus, L=0.0):
        Circuit.__init__(self, plus, minus, L=L)
        self.branches.append(Branch(plus, minus))
    def G(self, x, epar=defaultepar):
        return self._G
    def C(self, x, epar=defaultepar):
        n = self.n
        C = zeros((n,n), dtype=object)
        C[-1,-1] = -self.ipar.L
        return C

class VS(Circuit):
    """Independent voltage source

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> DC(c).solve(refnode=gnd).x
    array([ 1.5   ,  0.    , -0.0015])
    
    """
    terminals = ['plus', 'minus']
    instparams = [Parameter(name='v', desc='Source voltage', 
                            unit='V', default=1),
                  Parameter(name='vac', desc='AC analysis amplitude', 
                            unit='V', default=1),
                  Parameter(name='noisePSD', 
                            desc='Voltage noise power spectral density', 
                            unit='V^2/Hz', default=0)]

    def __init__(self, plus, minus, **kvargs):
        Circuit.__init__(self, plus, minus, **kvargs)
        self.branches.append(Branch(self.nodenames['plus'], 
                                    self.nodenames['minus']))

    def G(self, x, epar=defaultepar):
        return array([[0 , 0, 1],
                      [0 , 0, -1],
                      [1 , -1, 0]], dtype=object)

    def u(self, t=0.0, epar=defaultepar, analysis=None):
        if analysis == 'ac':
            return array([0, 0, -self.ipar.vac], dtype=object)
        elif analysis == None:
            return array([0, 0, -self.ipar.v], dtype=object)
        else:
            return super(VS, self).u(t,epar,analysis)

    def CY(self, x, w, epar=defaultepar):
        CY = super(VS, self).CY(x, w)
        CY[2, 2] = self.ipar.noisePSD
        return CY

    @property
    def branch(self):
        """Return the branch (plus, minus)"""
        return self.branches[0]

class IS(Circuit):
    """Independent DC current source

    >>> from analysis import DC, gnd as gnd2
    >>> c = SubCircuit()
    >>> n1=c.add_node('1')
    >>> c['is'] = IS(gnd, n1, i=1e-3)
    >>> c['R'] = R(n1, gnd, r=1e3)
    >>> DC(c).solve(refnode=gnd).x
    array([ 1.,  0.])

    """
    instparams = [Parameter(name='i', desc='DC Current', 
                            unit='A', default=1e-3),
                  Parameter(name='iac', desc='Small signal current amplitude', 
                            unit='A', default=0),
                  Parameter(name='noisePSD', 
                            desc='Current noise power spectral density', 
                            unit='A^2/Hz', default=0.0)]
    terminals = ['plus', 'minus']

    def u(self, t=0.0, epar=defaultepar, analysis=None):
        if analysis == None:
            return array([self.ipar.i, -self.ipar.i])
        elif analysis == 'ac':
            return array([self.ipar.iac, -self.ipar.iac])
        else:
            return super(IS, self).u(t,epar,analysis)

    def CY(self, x, w, epar=defaultepar):
        return  array([[self.ipar.noisePSD, -self.ipar.noisePSD],
                       [-self.ipar.noisePSD, self.ipar.noisePSD]], dtype=object)

class VCVS(Circuit):
    """Voltage controlled voltage source

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1, n2 =c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = VCVS(n1, gnd, n2, gnd, g=2)
    >>> c.nodes
    [Node('1'), Node('2'), Node('gnd', isglobal=True)]
    >>> c.branches
    [Branch(Node('1'),Node('gnd', isglobal=True)), Branch(Node('2'),Node('gnd', isglobal=True))]
    >>> c['vcvs'].G(zeros(4))
    array([[0, 0, 0, 0, 0],
           [0, 0, 0, 0, 0],
           [0, 0, 0, 0, 1],
           [0, 0, 0, 0, -1],
           [2, -2, -1, 1, 0]], dtype=object)

    """
    instparams = [Parameter(name='g', desc='Voltage gain', 
                            unit='V/V', default=1)]
    terminals = ['inp', 'inn', 'outp', 'outn']
    def __init__(self, *args, **kvargs):
        Circuit.__init__(self, *args, **kvargs)
        self.branches.append(Branch(self.nodenames['outp'], 
                                    self.nodenames['outn']))
        
    def G(self, x, epar=defaultepar):
        G = super(VCVS, self).G(x)
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
           (self.nodes.index(self.nodenames[name]) 
            for name in ('inp', 'inn', 'outp', 'outn'))
        G[outpindex, branchindex] += 1
        G[outnindex, branchindex] += -1
        G[branchindex, outpindex] += -1
        G[branchindex, outnindex] += 1
        G[branchindex, inpindex] += self.ipar.g
        G[branchindex, innindex] += -self.ipar.g
        return G

class VCCS(Circuit):
    """Voltage controlled current source

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1,n2 = c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vccs'] = VCCS(n1, gnd, n2, gnd, gm=1e-3)
    >>> c['rl'] = R(n2, gnd, r=1e3)
    >>> DC(c).solve(refnode=gnd).x
    array([ 1.5, -1.5,  0. ,  0. ])

    """
    terminals = ['inp', 'inn', 'outp', 'outn']
    instparams = [Parameter(name='gm', desc='Transconductance', 
                            unit='A/V', default=1e-3)]
    
    def G(self, x, epar=defaultepar):
        G = super(VCCS, self).G(x)
        gm=self.ipar.gm
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))
        G[outpindex, inpindex] += gm
        G[outpindex, innindex] += -gm
        G[outnindex, inpindex] += -gm
        G[outnindex, innindex] += gm
        return G

class Nullor(Circuit):
    """Nullor

    From Wikipedia, the free encyclopedia

     A nullor is a theoretical two-port network comprised of a nullator at 
     its input and a norator at its output.[1]
     Nullors represent an ideal amplifier, having infinite current, 
     voltage, transconductance and transimpedance gain.[2] 
     Its transmission parameters are all zero.

     1. The name "nullor" was introduced by H.J. Carlin, 
      Singular network elements, 
      IEEE Trans. Circuit Theory, March 1965, vol. CT-11, pp. 67-72.
 
     2. Verhoeven C J M van Staveren A Monna G L E Kouwenhoven, 
       M H L & Yildiz E (2003). 
       Structured electronic design: negative feedback amplifiers.
       Boston/Dordrecht/London: Kluwer Academic, §2.2.2 pp. 32-34. 
       ISBN 1402075901.

    """
    terminals = ('inp', 'inn', 'outp', 'outn')

    def __init__(self, *args, **kvargs):
        Circuit.__init__(self, *args, **kvargs)
        self.branches.append(Branch(self.nodenames['outp'], 
                                    self.nodenames['outn']))

    def G(self, x, epar=defaultepar):
        G = super(Nullor, self).G(x)
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))

        G[outpindex, branchindex] += 1
        G[outnindex, branchindex] += -1
        G[branchindex, inpindex] += 1
        G[branchindex, innindex] += -1
        return G


class Transformer(Circuit):
    """Ideal transformer

    >>> from analysis import DC
    >>> c = SubCircuit()
    >>> n1, n2 = c.add_nodes('1', '2')
    >>> c['vs'] = VS(n1, gnd, v=1.5)
    >>> c['vcvs'] = Transformer(n1, gnd, n2, gnd, n=2)
    >>> c['vcvs'].nodes
    [Node('inp'), Node('inn'), Node('outp'), Node('outn')]
    >>> c['vcvs'].branches
    [Branch(Node('outp'),Node('outn'))]
    >>> c['vcvs'].G(zeros(4))
    array([[0, 0, 0, 0, 2],
           [0, 0, 0, 0, -2],
           [0, 0, 0, 0, 1],
           [0, 0, 0, 0, -1],
           [-1, 1, 2, -2, 0]], dtype=object)

    """
    instparams = [Parameter(name='n', desc='Winding ratio', unit='', default=1)]
    terminals = ['inp', 'inn', 'outp', 'outn']
    def __init__(self, *args, **kvargs):
        Circuit.__init__(self, *args, **kvargs)
        self.branches.append(Branch(self.nodenames['outp'], 
                                    self.nodenames['outn']))

    def G(self, x, epar=defaultepar):
        G = super(Transformer, self).G(x)
        branchindex = -1
        inpindex, innindex, outpindex, outnindex = \
            (self.nodes.index(self.nodenames[name]) 
             for name in ('inp', 'inn', 'outp', 'outn'))
        G[inpindex, branchindex] += self.ipar.n
        G[innindex, branchindex] += -self.ipar.n
        G[outpindex, branchindex] += 1
        G[outnindex, branchindex] += -1
        G[branchindex, outpindex] += self.ipar.n
        G[branchindex, outnindex] += -self.ipar.n
        G[branchindex, inpindex] += -1
        G[branchindex, innindex] += 1
        return G

class Diode(Circuit):
    terminals = ['plus', 'minus']
    mpar = Circuit.mpar.copy( 
        Parameter(name='IS', desc='Saturation current', 
                  unit='A', default=1e-13))
    linear = False
    def G(self, x, epar=defaultepar):
        VD = x[0]-x[1]
        VT = kboltzmann*epar.T / qelectron
        g = self.mpar.IS*exp(VD/VT)/VT
        return array([[g, -g],
                      [-g, g]], dtype=object)

    def i(self, x, epar=defaultepar):
        """
        
        """
        VD = x[0]-x[1]
        VT = kboltzmann*epar.T / qelectron
        I = self.mpar.IS*(exp(VD/VT)-1.0)
        return array([I, -I])



if __name__ == "__main__":
    import doctest
    doctest.testmod()
